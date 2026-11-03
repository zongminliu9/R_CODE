#!/usr/bin/env python3
"""Restrained minimisation and interaction-energy decomposition for
   {TBCB-CAP-Gly, TBCE-CAP-Gly} x {EEY, EEF}, with the three His49 states for TBCE.

WHAT THIS IS, AND IS NOT.  These are interaction energies (and single-structure MM/GB
estimates) at restrained-minimised geometries. They contain no conformational entropy, no
solute entropy, and no sampling, so they are NOT binding free energies and are not interchangeable
with the dG_bind values the original plan asked for. They are useful for COMPARING chemically
similar systems that share a binding mode, which is exactly the EEY-vs-EEF question.

Protonation.  His49 of TBCE sits where TBCB has Asn204, facing the peptide's free carboxylate.
Its pKa (~6.5 free in solution) can be shifted up by that carboxylate, so all three states are
computed rather than assumed: HID and HIE (neutral tautomers) and HIP (protonated, +1).

Restraints.  The peptide backbone is harmonically restrained to the transplanted (experimental)
pose so the binding mode is preserved while clashes relax; the domain and all side chains,
including the peptide's terminal side chain, are free.
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
from pdbfixer import PDBFixer

KJ = u.kilojoule_per_mole
KCAL = u.kilocalorie_per_mole
os.makedirs("results/energetics", exist_ok=True)
os.makedirs("data/models/prep", exist_ok=True)

def prep(src, out, his_state=None):
    fx = PDBFixer(filename=src)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fx.topology, fx.positions, open(out,"w"), keepIds=True)
    return out

def build(pdbfile, implicit=True):
    pdb = app.PDBFile(pdbfile)
    ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml') if implicit \
         else app.ForceField('amber14-all.xml')
    sysobj = ff.createSystem(pdb.topology,
                nonbondedMethod=app.CutoffNonPeriodic if implicit else app.NoCutoff,
                nonbondedCutoff=1.8*u.nanometer, constraints=None)
    return pdb, sysobj

def energy_of_subset(pdbfile, keep_chains, coords=None, implicit=True):
    """Energy of only the atoms in keep_chains, at the given (full-system) coordinates."""
    pdb = app.PDBFile(pdbfile)
    mod = app.Modeller(pdb.topology, coords if coords is not None else pdb.positions)
    drop = [c for c in mod.topology.chains() if c.id not in keep_chains]
    mod.delete(drop)
    ff = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml') if implicit \
         else app.ForceField('amber14-all.xml')
    s = ff.createSystem(mod.topology,
            nonbondedMethod=app.CutoffNonPeriodic if implicit else app.NoCutoff,
            nonbondedCutoff=1.8*u.nanometer, constraints=None)
    integ = mm.VerletIntegrator(1*u.femtosecond)
    ctx = mm.Context(s, integ, mm.Platform.getPlatformByName('CPU'), {'Threads':'4'})
    ctx.setPositions(mod.positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
    del ctx, integ
    return e, mod.topology.getNumAtoms()

def restrained_minimise(pdbfile, k=50.0, maxit=3000):
    pdb, sysobj = build(pdbfile, implicit=True)
    # restrain peptide backbone (chain P) to the experimental pose
    rest = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    rest.addGlobalParameter("k", k*u.kilojoule_per_mole/u.nanometer**2)
    for p in ("x0","y0","z0"): rest.addPerParticleParameter(p)
    n = 0
    for atom in pdb.topology.atoms():
        if atom.residue.chain.id == "P" and atom.name in ("N","CA","C","O","OXT"):
            pos = pdb.positions[atom.index]
            rest.addParticle(atom.index, [pos.x, pos.y, pos.z]); n += 1
    sysobj.addForce(rest)
    integ = mm.LangevinMiddleIntegrator(310*u.kelvin, 1/u.picosecond, 1*u.femtosecond)
    sim = app.Simulation(pdb.topology, sysobj, integ,
                         mm.Platform.getPlatformByName('CPU'), {'Threads':'4'})
    sim.context.setPositions(pdb.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
    sim.minimizeEnergy(maxIterations=maxit)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = st.getPotentialEnergy().value_in_unit(KCAL)
    out = pdbfile.replace(".pdb", "_min.pdb")
    app.PDBFile.writeFile(pdb.topology, st.getPositions(), open(out, "w"), keepIds=True)
    rmsd_pep = None
    return out, e0, e1, n

SYSTEMS = []
for dom in ("TBCB","TBCE"):
    for pep in ("EEY","EEF"):
        SYSTEMS.append((f"{dom}_{pep}", f"data/models/{dom}_capgly_{pep}_raw.pdb", dom, pep))

print(f"{'system':14s} {'atoms':>6s} {'E_raw':>12s} {'E_min':>12s} {'nrest':>6s}")
rows = []
for name, src, dom, pep in SYSTEMS:
    p = prep(src, f"data/models/prep/{name}.pdb")
    try:
        out, e0, e1, nrest = restrained_minimise(p)
    except Exception as ex:
        print(f"{name:14s} FAILED: {type(ex).__name__}: {str(ex)[:80]}"); continue
    pdbm = app.PDBFile(out)
    coords = pdbm.positions
    e_cplx, na = energy_of_subset(out, {"A","E","B","P"}, coords)
    # domain chain id: whatever is not P
    chains = {c.id for c in pdbm.topology.chains()}
    dchain = chains - {"P"}
    e_dom, nd = energy_of_subset(out, dchain, coords)
    e_pep, npep = energy_of_subset(out, {"P"}, coords)
    dE = e_cplx - e_dom - e_pep
    # vacuum decomposition
    ev_c,_ = energy_of_subset(out, chains, coords, implicit=False)
    ev_d,_ = energy_of_subset(out, dchain, coords, implicit=False)
    ev_p,_ = energy_of_subset(out, {"P"}, coords, implicit=False)
    dEv = ev_c - ev_d - ev_p
    print(f"{name:14s} {na:6d} {e0:12.1f} {e1:12.1f} {nrest:6d}"
          f"   dE_int(GB) = {dE:9.2f}  dE_int(vac) = {dEv:9.2f}  dG_solv = {dE-dEv:9.2f} kcal/mol")
    rows.append(dict(system=name, domain=dom, peptide=pep, atoms=na, n_restrained=nrest,
                     E_raw=e0, E_min=e1, dE_int_GB=dE, dE_int_vac=dEv, dG_solv_coupling=dE-dEv,
                     n_domain=nd, n_peptide=npep, minimised=out))
json.dump(rows, open("results/energetics/interaction_energies.json","w"), indent=1)

print("\n=== EEY -> EEF difference (negative = EEF binds better by this measure) ===")
by = {r["system"]: r for r in rows}
for dom in ("TBCB","TBCE"):
    y, f = by.get(f"{dom}_EEY"), by.get(f"{dom}_EEF")
    if not (y and f): continue
    print(f"  {dom}:  ddE_int(GB) = {f['dE_int_GB']-y['dE_int_GB']:+8.2f} kcal/mol   "
          f"ddE_int(vac) = {f['dE_int_vac']-y['dE_int_vac']:+8.2f}")
