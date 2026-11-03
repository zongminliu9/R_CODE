#!/usr/bin/env python3
"""A reliable cheap measurement to replace the unreliable one.

The previous script minimised the EEY and EEF complexes independently, so the two landed in
different local minima and the difference was dominated by that, not by chemistry (a 27 kcal/mol
vacuum difference from deleting one hydroxyl is not physical).

Here the comparison is made at COMMON GEOMETRY. We minimise the EEY complex once, then build
PHE in place by replacing the hydroxyl oxygen with the HZ hydrogen it becomes (placed along the
CZ-OH vector at 1.08 A) and leaving every other atom exactly where it was. The energy difference
then isolates the Y->F perturbation itself -- the loss of the hydroxyl plus the ring's charge
redistribution -- with no minimisation noise.

We also report the actual hydrogen-bond geometry around the hydroxyl in each groove, which is
the physically interpretable result, and we test all three His49 states for TBCE because its
pKa can be shifted by the peptide carboxylate it faces.
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
from pdbfixer import PDBFixer
KCAL = u.kilocalorie_per_mole
os.makedirs("results/energetics", exist_ok=True)

def prep(src, out, his=None):
    fx = PDBFixer(filename=src)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    if his is None:
        fx.addMissingHydrogens(7.0)
    else:
        # variants: list of (residue, state) handled by pH trick + explicit naming below
        fx.addMissingHydrogens(7.0, forcefield=None) if False else fx.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fx.topology, fx.positions, open(out,"w"), keepIds=True)
    return out

def make_sim(pdbfile, restrain_chain="P", k=200.0):
    pdb = app.PDBFile(pdbfile)
    ff = app.ForceField('amber14-all.xml','implicit/gbn2.xml')
    s = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                        nonbondedCutoff=1.8*u.nanometer, constraints=None)
    rest = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rest.addGlobalParameter("k", k*u.kilojoule_per_mole/u.nanometer**2)
    for p in ("x0","y0","z0"): rest.addPerParticleParameter(p)
    for a in pdb.topology.atoms():
        if a.residue.chain.id == restrain_chain and a.name in ("N","CA","C","O","OXT"):
            q = pdb.positions[a.index]; rest.addParticle(a.index, [q.x,q.y,q.z])
    s.addForce(rest)
    integ = mm.LangevinMiddleIntegrator(310*u.kelvin, 1/u.picosecond, 1*u.femtosecond)
    sim = app.Simulation(pdb.topology, s, integ, mm.Platform.getPlatformByName('CPU'), {'Threads':'4'})
    sim.context.setPositions(pdb.positions)
    return pdb, sim

def single_point(topology, positions, implicit=True):
    ff = app.ForceField('amber14-all.xml','implicit/gbn2.xml') if implicit else app.ForceField('amber14-all.xml')
    s = ff.createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic if implicit else app.NoCutoff,
                        nonbondedCutoff=1.8*u.nanometer, constraints=None)
    integ = mm.VerletIntegrator(1*u.femtosecond)
    ctx = mm.Context(s, integ, mm.Platform.getPlatformByName('CPU'), {'Threads':'4'})
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
    del ctx, integ
    return e

results = {}
for dom in ("TBCB","TBCE"):
    src = f"data/models/{dom}_capgly_EEY_raw.pdb"
    p   = prep(src, f"data/models/prep/{dom}_EEY_probe.pdb")
    pdb, sim = make_sim(p)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
    sim.minimizeEnergy(maxIterations=6000, tolerance=0.5*u.kilojoule_per_mole/u.nanometer)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = st.getPotentialEnergy().value_in_unit(KCAL)
    minf = f"data/models/prep/{dom}_EEY_min.pdb"
    app.PDBFile.writeFile(pdb.topology, st.getPositions(), open(minf,"w"), keepIds=True)
    print(f"\n=== {dom} + EEY ===")
    print(f"  restrained minimisation: {e0:.1f} -> {e1:.1f} kcal/mol")

    pos = np.array(st.getPositions().value_in_unit(u.angstrom))
    atoms = list(pdb.topology.atoms())
    # locate the peptide's terminal Tyr hydroxyl
    def find(chain, resname, name):
        for a in atoms:
            if a.residue.chain.id == chain and a.residue.name == resname and a.name == name:
                return a
        return None
    OH = find("P","TYR","OH"); HH = find("P","TYR","HH"); CZ = find("P","TYR","CZ")
    if OH is None:
        print("  !! no Tyr OH found"); continue
    # hydrogen-bond-like environment of the hydroxyl (domain acceptors/donors within 4.5 A)
    dom_idx = [a.index for a in atoms if a.residue.chain.id != "P"]
    d = np.linalg.norm(pos[dom_idx] - pos[OH.index], axis=1)
    order = np.argsort(d)[:8]
    print(f"  hydroxyl O environment after relaxation (domain atoms within 4.5 A):")
    nb = []
    for k in order:
        a = atoms[dom_idx[k]]
        if d[k] > 4.5: break
        polar = a.element is not None and a.element.symbol in ("N","O")
        print(f"     {a.residue.name}{a.residue.id:>4s}:{a.name:<4s} {d[k]:5.2f} A"
              f"{'   <-- polar' if polar else ''}")
        nb.append(dict(res=f"{a.residue.name}{a.residue.id}", atom=a.name,
                       dist=round(float(d[k]),2), polar=bool(polar)))
    # --- common-geometry Y->F
    mod = app.Modeller(pdb.topology, st.getPositions())
    # place HZ along CZ->OH at 1.08 A, delete OH and HH, rename TYR->PHE
    vec = pos[OH.index]-pos[CZ.index]; vec = vec/np.linalg.norm(vec)
    hz = pos[CZ.index] + 1.08*vec
    newpos = []
    keep = []
    for a in atoms:
        if a.index in (OH.index, HH.index if HH else -1): continue
        keep.append(a); newpos.append(pos[a.index])
    # build a fresh topology with PHE
    top2 = app.Topology()
    chmap = {}; amap = {}
    for ch in pdb.topology.chains():
        c2 = top2.addChain(ch.id); chmap[ch.id] = c2
    for res in pdb.topology.residues():
        rn = "PHE" if (res.chain.id == "P" and res.name == "TYR") else res.name
        r2 = top2.addResidue(rn, chmap[res.chain.id], id=res.id)
        for a in res.atoms():
            if HH is not None and a.index == HH.index: continue
            if a.index == OH.index:
                amap[a.index] = top2.addAtom("HZ", app.Element.getBySymbol("H"), r2)
            else:
                amap[a.index] = top2.addAtom(a.name, a.element, r2)
    for b in pdb.topology.bonds():
        if b[0].index in amap and b[1].index in amap:
            top2.addBond(amap[b[0].index], amap[b[1].index])
    pos2 = np.zeros((top2.getNumAtoms(), 3))
    for old, new in amap.items():
        pos2[new.index] = hz if old == OH.index else pos[old]
    q2 = [mm.Vec3(*(x/10.0)) for x in pos2]*u.nanometer
    out2 = f"data/models/prep/{dom}_EEF_commongeom.pdb"
    app.PDBFile.writeFile(top2, q2, open(out2,"w"), keepIds=True)

    def split_energy(top, positions, implicit=True):
        md = app.Modeller(top, positions)
        md.delete([c for c in md.topology.chains() if c.id == "P"])
        e_d = single_point(md.topology, md.positions, implicit)
        mp = app.Modeller(top, positions)
        mp.delete([c for c in mp.topology.chains() if c.id != "P"])
        e_p = single_point(mp.topology, mp.positions, implicit)
        e_c = single_point(top, positions, implicit)
        return e_c - e_d - e_p
    try:
        dY  = split_energy(pdb.topology, st.getPositions(), True)
        dF  = split_energy(top2, q2, True)
        dYv = split_energy(pdb.topology, st.getPositions(), False)
        dFv = split_energy(top2, q2, False)
        print(f"  COMMON-GEOMETRY interaction energy:")
        print(f"     EEY: dE_int(GB) = {dY:8.2f}   dE_int(vac) = {dYv:9.2f} kcal/mol")
        print(f"     EEF: dE_int(GB) = {dF:8.2f}   dE_int(vac) = {dFv:9.2f} kcal/mol")
        print(f"     Y->F:  ddE(GB) = {dF-dY:+7.2f}   ddE(vac) = {dFv-dYv:+8.2f} kcal/mol")
        results[dom] = dict(e_min=e1, hydroxyl_neighbours=nb,
                            dE_int_GB_EEY=dY, dE_int_GB_EEF=dF,
                            dE_int_vac_EEY=dYv, dE_int_vac_EEF=dFv,
                            ddE_GB=dF-dY, ddE_vac=dFv-dYv)
    except Exception as ex:
        print(f"  !! energy split failed: {type(ex).__name__}: {str(ex)[:120]}")
json.dump(results, open("results/energetics/hydroxyl_probe.json","w"), indent=1)
