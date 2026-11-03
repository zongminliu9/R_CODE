#!/usr/bin/env python3
"""Build, solvate and equilibrate every system the GPU jobs need.

Deliberate scope note. This script does the parts that can be written and checked without a
GPU: structure preparation, protonation, solvation, parametrisation, minimisation and a short
restrained equilibration. The alchemical hybrid-topology construction is NOT hand-rolled here --
10_fep_mutation.py delegates it to perses, because mapping Ile->Ser correctly (CB-CG1-CG2-CD1
against CB-OG) is the step where a home-made implementation silently produces a confident wrong
number. Everything this file writes is a plain OpenMM System/State, so it is inspectable.

Outputs, per system, into systems/:
    <name>_solv.pdb     equilibrated coordinates
    <name>_solv.xml     serialised System
    <name>_state.xml    serialised State (positions, velocities, box)
    <name>_meta.json    composition, charge, atom counts, build settings
"""
import argparse, json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
from pdbfixer import PDBFixer

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT  = os.path.join(HERE, "systems")

TEMP      = 310*u.kelvin
SALT      = 0.15*u.molar
PADDING   = 1.2*u.nanometer
FF_PROT   = "amber14-all.xml"     # ff19SB where available: see --ff19sb
FF_WATER  = "amber14/opc.xml"
HMASS     = 3.0*u.amu
DT_EQ     = 2.0*u.femtoseconds

def forcefield(ff19sb: bool):
    """ff19SB + OPC is what the study design calls for. OpenMM ships it in recent
    openmmforcefields releases; fall back to ff14SB + OPC and SAY SO rather than
    silently changing the Hamiltonian."""
    if ff19sb:
        for name in ("amber19-all.xml", "amber/ff19SB.xml", "amber19/ff19SB.xml"):
            try:
                ff = app.ForceField(name, FF_WATER); return ff, name
            except Exception:
                continue
        print("  !! ff19SB not found in this OpenMM install; falling back to "
              "amber14-all.xml (ff14SB). Record this in the methods.", file=sys.stderr)
    return app.ForceField(FF_PROT, FF_WATER), FF_PROT

def prep_structure(src, out_pdb, keep_chains=None, res_range=None, ph=7.0):
    fx = PDBFixer(filename=src)
    if keep_chains is not None:
        drop = [i for i, c in enumerate(fx.topology.chains()) if c.id not in keep_chains]
        if drop: fx.removeChains(drop)
    fx.findMissingResidues()
    # do not build unresolved terminal residues: they would be modelled, not measured
    fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(ph)
    app.PDBFile.writeFile(fx.topology, fx.positions, open(out_pdb, "w"), keepIds=True)
    return out_pdb

def build(name, pdb_in, ff19sb=False, equil_ps=200.0, minimise_only=False):
    os.makedirs(OUT, exist_ok=True)
    ff, ffname = forcefield(ff19sb)
    pdb = app.PDBFile(pdb_in)
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.addSolvent(ff, model="tip4pew", padding=PADDING, ionicStrength=SALT, neutralize=True)
    n_atoms = mod.topology.getNumAtoms()
    sysobj = ff.createSystem(mod.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0*u.nanometer, constraints=app.HBonds,
                             rigidWater=True, hydrogenMass=HMASS)
    sysobj.addForce(mm.MonteCarloBarostat(1*u.bar, TEMP, 25))
    integ = mm.LangevinMiddleIntegrator(TEMP, 1/u.picosecond, DT_EQ)
    plat = mm.Platform.getPlatformByName(os.environ.get("OPENMM_DEFAULT_PLATFORM", "CPU"))
    sim = app.Simulation(mod.topology, sysobj, integ, plat)
    sim.context.setPositions(mod.positions)
    # Measure the starting energy AFTER applying constraints. minimizeEnergy enforces the
    # HBonds and rigid-water constraints, so an unconstrained e0 is not comparable with the
    # constrained e1 and makes minimisation look as though it raised the energy.
    sim.context.applyConstraints(1e-8)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=5000)
    e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  {name}: {n_atoms} atoms, ff={ffname}; minimised "
          f"{e0.value_in_unit(u.kilocalorie_per_mole):.1f} -> "
          f"{e1.value_in_unit(u.kilocalorie_per_mole):.1f} kcal/mol", flush=True)
    if not minimise_only:
        # restrained heating on the solute, then free equilibration
        rest = mm.CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
        rest.addGlobalParameter("k", 10.0*u.kilocalorie_per_mole/u.angstrom**2)
        for p in ("x0","y0","z0"): rest.addPerParticleParameter(p)
        solute = [a.index for a in mod.topology.atoms()
                  if a.residue.name not in ("HOH","NA","CL") and a.element.symbol != "H"]
        pos = sim.context.getState(getPositions=True).getPositions()
        for i in solute:
            rest.addParticle(i, [pos[i].x, pos[i].y, pos[i].z])
        fid = sysobj.addForce(rest); sim.context.reinitialize(preserveState=True)
        sim.context.setVelocitiesToTemperature(TEMP)
        sim.step(int(equil_ps*0.25*1000/DT_EQ.value_in_unit(u.femtoseconds)))
        sysobj.removeForce(fid); sim.context.reinitialize(preserveState=True)
        sim.step(int(equil_ps*0.75*1000/DT_EQ.value_in_unit(u.femtoseconds)))
    st = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    app.PDBFile.writeFile(mod.topology, st.getPositions(), open(f"{OUT}/{name}_solv.pdb","w"))
    open(f"{OUT}/{name}_solv.xml","w").write(mm.XmlSerializer.serialize(sysobj))
    open(f"{OUT}/{name}_state.xml","w").write(mm.XmlSerializer.serialize(st))
    charge = sum(sysobj.getForce(i).getParticleParameters(j)[0].value_in_unit(u.elementary_charge)
                 for i in range(sysobj.getNumForces())
                 if isinstance(sysobj.getForce(i), mm.NonbondedForce)
                 for j in range(sysobj.getForce(i).getNumParticles()))
    meta = dict(name=name, atoms=n_atoms, forcefield=ffname, water="OPC/tip4pew",
                temperature_K=TEMP.value_in_unit(u.kelvin),
                ionic_strength_M=SALT.value_in_unit(u.molar),
                padding_nm=PADDING.value_in_unit(u.nanometer),
                hydrogen_mass_amu=HMASS.value_in_unit(u.amu),
                equil_ps=0.0 if minimise_only else equil_ps,
                net_charge_e=round(float(charge), 4),
                source_pdb=os.path.relpath(pdb_in, HERE))
    json.dump(meta, open(f"{OUT}/{name}_meta.json","w"), indent=1)
    return meta

def build_fragment(name, src_pdb, lo, hi, ff19sb=False):
    """Unfolded reference for a ddG_fold cycle: the mutation site in its own local sequence,
    excised from a real structure so the geometry is physical. Hand-building an extended
    conformation by placing atoms along an axis gives colinear overlaps and NaNs -- that bug
    was hit and is why this reads coordinates instead.

    Residues lo..hi are taken with the site in the middle; standard charged termini are used,
    and because the same fragment appears in both legs of the cycle their contribution cancels.
    """
    import biotite.structure as struc, biotite.structure.io.pdb as bpdb
    f = bpdb.PDBFile.read(src_pdb)
    arr = f.get_structure(model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    sub = arr[(arr.res_id >= lo) & (arr.res_id <= hi)]
    if sub.array_length() == 0:
        raise ValueError(f"no residues {lo}-{hi} in {src_pdb}")
    seq = "".join(sorted({(int(r), n) for r, n in zip(sub.res_id, sub.res_name)}) and
                  [n for _, n in sorted({(int(r), n) for r, n in zip(sub.res_id, sub.res_name)})])
    raw = f"{OUT}/_{name}_raw.pdb"
    out = bpdb.PDBFile(); out.set_structure(sub); out.write(raw)
    print(f"  fragment {lo}-{hi} from {os.path.basename(src_pdb)}: "
          f"{len(set(sub.res_id.tolist()))} residues, {sub.array_length()} atoms")
    prep_structure(raw, f"{OUT}/_{name}_prep.pdb")
    return build(name, f"{OUT}/_{name}_prep.pdb", ff19sb=ff19sb, equil_ps=100.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--ff19sb", action="store_true", default=True)
    ap.add_argument("--minimise-only", action="store_true",
                    help="skip equilibration; for checking the build path on CPU")
    ap.add_argument("--equil-ps", type=float, default=200.0)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    AF3 = "/home/user/uploads_fold"
    jobs = {
      # unfolded reference for I114S: residues 112-116 of TUBA1A, site in the middle
      "fragment_I114": lambda: build_fragment("fragment_I114",
          f"{REPO}/data/afdb/P68369.pdb", 112, 116, a.ff19sb),
      "tuba1a": lambda: build("tuba1a",
          prep_structure(f"{REPO}/data/afdb/P68369.pdb", f"{OUT}/_tuba1a_prep.pdb"),
          a.ff19sb, a.equil_ps, a.minimise_only),
      "tuba8": lambda: build("tuba8",
          prep_structure(f"{REPO}/data/afdb/Q9JJZ2.pdb", f"{OUT}/_tuba8_prep.pdb"),
          a.ff19sb, a.equil_ps, a.minimise_only),
      "tbcb_tuba1a": lambda: build("tbcb_tuba1a",
          prep_structure(f"{AF3}/88398339-fold_tbcb_tuba1a/fold_tbcb_tuba1a_model_0.cif",
                         f"{OUT}/_tbcb_tuba1a_prep.pdb"),
          a.ff19sb, a.equil_ps, a.minimise_only),
      "tbcb_tuba8": lambda: build("tbcb_tuba8",
          prep_structure(f"{AF3}/e50dbe95-fold_tbcb_tuba8/fold_tbcb_tuba8_model_0.cif",
                         f"{OUT}/_tbcb_tuba8_prep.pdb"),
          a.ff19sb, a.equil_ps, a.minimise_only),
    }
    want = list(jobs) if (a.all or not a.only) else a.only
    built = {}
    for k in want:
        if k not in jobs: print(f"  unknown system {k}", file=sys.stderr); continue
        print(f"building {k} ...", flush=True)
        try:
            built[k] = jobs[k]()
        except Exception as e:
            print(f"  !! {k} failed: {type(e).__name__}: {e}", file=sys.stderr)
    json.dump(built, open(f"{OUT}/manifest.json","w"), indent=1)
    print(f"\n{len(built)}/{len(want)} systems built -> {OUT}")

if __name__ == "__main__":
    main()
