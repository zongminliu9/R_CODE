#!/usr/bin/env python3
"""Benchmark OpenMM CPU throughput on the ACTUAL systems this project needs."""
import time, sys, numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
from pdbfixer import PDBFixer

def prep(cif, keep_chains, out, ph=7.0):
    fx = PDBFixer(filename=cif)
    # drop everything except requested author chains
    allch = [c.id for c in fx.topology.chains()]
    drop = [i for i,c in enumerate(fx.topology.chains()) if c.id not in keep_chains]
    fx.removeChains(drop)
    fx.findMissingResidues()
    # do not build un-resolved terminal/loop residues at chain ends for a benchmark
    fx.missingResidues = {}
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms(); fx.addMissingAtoms()
    fx.addMissingHydrogens(ph)
    app.PDBFile.writeFile(fx.topology, fx.positions, open(out,"w"), keepIds=True)
    n = sum(1 for _ in fx.topology.atoms()); nres = sum(1 for _ in fx.topology.residues())
    return out, n, nres

def bench(pdb_path, label, explicit, nsteps=800, dt=0.004):
    pdb = app.PDBFile(pdb_path)
    if explicit:
        ff = app.ForceField('amber14-all.xml','amber14/opc.xml')
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.addSolvent(ff, model='tip4pew', padding=1.2*u.nanometer,
                       ionicStrength=0.15*u.molar, neutralize=True)
        top, pos = mod.topology, mod.positions
        sysobj = ff.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0*u.nanometer,
                                 constraints=app.HBonds, rigidWater=True, hydrogenMass=None)
        sysobj.addForce(mm.MonteCarloBarostat(1*u.bar, 310*u.kelvin, 50))
    else:
        ff = app.ForceField('amber14-all.xml','implicit/gbn2.xml')
        top, pos = pdb.topology, pdb.positions
        sysobj = ff.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=2.0*u.nanometer, constraints=app.HBonds,
                                 soluteDielectric=1.0, solventDielectric=78.5)
    integ = mm.LangevinMiddleIntegrator(310*u.kelvin, 1/u.picosecond, dt*u.picoseconds)
    sim = app.Simulation(top, sysobj, integ, mm.Platform.getPlatformByName('CPU'), {'Threads':'4'})
    sim.context.setPositions(pos)
    sim.minimizeEnergy(maxIterations=400)
    sim.context.setVelocitiesToTemperature(310*u.kelvin)
    sim.step(200)
    t0 = time.time(); sim.step(nsteps); el = time.time()-t0
    nsday = (nsteps*dt/1000.0)/el*86400
    na = sysobj.getNumParticles()
    print(f"{label:42s} {na:7d} atoms  {nsday:8.2f} ns/day   100ns = {100/nsday:7.1f} d")
    sys.stdout.flush()
    return nsday, na

print("preparing real systems ...")
p1,n1,r1 = prep("data/pdb/2HL3.cif", {"A","C"}, "data/pdb/prep_2HL3_capgly_EEY.pdb")
print(f"  2HL3 CAP-Gly+EEY : {r1} res, {n1} atoms")
p2,n2,r2 = prep("data/pdb/1JFF.cif", {"A"}, "data/pdb/prep_1JFF_alpha.pdb")
print(f"  1JFF alpha-tubulin: {r2} res, {n2} atoms")
p3,n3,r3 = prep("data/pdb/1WHG.cif", {"A"}, "data/pdb/prep_1WHG_mTBCB_capgly.pdb")
print(f"  1WHG mTBCB CAP-Gly: {r3} res, {n3} atoms\n")

print(f"{'system':42s} {'atoms':>7s}         {'throughput':>8s}")
print("-"*92)
bench(p1, "CAP-Gly+EEY   (2HL3)  GBn2 implicit", False)
bench(p3, "mTBCB CAP-Gly (1WHG)  GBn2 implicit", False)
bench(p2, "alpha-tubulin (1JFF)  GBn2 implicit", False)
bench(p1, "CAP-Gly+EEY   (2HL3)  OPC explicit", True, nsteps=400)
bench(p2, "alpha-tubulin (1JFF)  OPC explicit", True, nsteps=200)
