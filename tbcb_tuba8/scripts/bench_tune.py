#!/usr/bin/env python3
"""Can cutoff / threading / timestep choices buy back an order of magnitude on CPU?"""
import time, sys
import openmm as mm, openmm.app as app, openmm.unit as u

def one(pdbf, cutoff_nm, threads, dt_fs, hmr, nsteps=600, label=""):
    pdb = app.PDBFile(pdbf)
    ff = app.ForceField('amber14-all.xml','implicit/gbn2.xml')
    kw = {}
    if hmr: kw['hydrogenMass'] = 3.0*u.amu
    if cutoff_nm is None:
        sysobj = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds, **kw)
        cl = "none"
    else:
        sysobj = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=cutoff_nm*u.nanometer,
                                 constraints=app.HBonds, **kw)
        cl = f"{cutoff_nm:.1f}nm"
    integ = mm.LangevinMiddleIntegrator(310*u.kelvin, 1/u.picosecond, dt_fs*0.001*u.picoseconds)
    sim = app.Simulation(pdb.topology, sysobj, integ,
                         mm.Platform.getPlatformByName('CPU'), {'Threads': str(threads)})
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=150)
    sim.context.setVelocitiesToTemperature(310*u.kelvin)
    sim.step(150)
    t0 = time.time(); sim.step(nsteps); el = time.time()-t0
    nsday = (nsteps*dt_fs*1e-6)/el*86400
    print(f"  {label:30s} cut={cl:6s} thr={threads} dt={dt_fs}fs hmr={int(hmr)}  "
          f"{nsday:7.2f} ns/day", flush=True)
    return nsday

P = "data/pdb/prep_2HL3_capgly_EEY.pdb"   # 1227 atoms -- the FEP-relevant system
print("=== CAP-Gly+EEY (1227 atoms): tuning ===", flush=True)
base = one(P, 2.0, 4, 4, False, label="baseline")
for cut in (1.5, 1.2, 1.0):
    one(P, cut, 4, 4, False, label=f"shorter GB cutoff")
one(P, 1.2, 4, 5, True,  label="cutoff1.2 + HMR dt=5fs")
print("\n=== threading scaling (for running lambda windows in parallel) ===", flush=True)
for thr in (1, 2, 4):
    one(P, 1.2, thr, 5, True, label=f"{thr}-thread")
