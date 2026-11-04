#!/usr/bin/env python3
"""Production MD of a CAP-Gly + C-terminal-peptide complex, with checkpointing.

One job = one system, one thread, so four run concurrently on the four cores (measured: 8.75
ns/day each, ~35 ns/day aggregate, against 26 ns/day for a single 4-thread job).

Checkpoints every 100 ps because this container has already restarted twice; a restart resumes
from the .chk rather than starting over. Frames every 10 ps for MM/GBSA ensemble averaging.
"""
import argparse, os, sys, time
import openmm as mm, openmm.app as app, openmm.unit as u

ap = argparse.ArgumentParser()
ap.add_argument("--pdb", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--ns", type=float, default=20.0)
ap.add_argument("--threads", default="1")
ap.add_argument("--implicit", default="implicit/gbn2.xml")
ap.add_argument("--cutoff", type=float, default=1.8)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

OUT = "md"; os.makedirs(OUT, exist_ok=True)
chk = f"{OUT}/{a.tag}.chk"; dcd = f"{OUT}/{a.tag}.dcd"; log = f"{OUT}/{a.tag}.csv"

pdb = app.PDBFile(a.pdb)
ff = app.ForceField("amber14-all.xml", a.implicit)
system = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                         nonbondedCutoff=a.cutoff*u.nanometer, constraints=app.HBonds,
                         hydrogenMass=3.0*u.amu, soluteDielectric=1.0, solventDielectric=78.5)
integ = mm.LangevinMiddleIntegrator(310*u.kelvin, 1/u.picosecond, 0.004*u.picoseconds)
integ.setRandomNumberSeed(a.seed)
sim = app.Simulation(pdb.topology, system, integ,
                     mm.Platform.getPlatformByName("CPU"), {"Threads": a.threads})
dt_ps = 0.004
steps_total = int(a.ns*1000/dt_ps)
resume = 0
if os.path.exists(chk):
    try:
        sim.loadCheckpoint(chk)
        resume = sim.currentStep
        print(f"[{a.tag}] resumed at step {resume} ({resume*dt_ps/1000:.3f} ns)", flush=True)
    except Exception as e:
        print(f"[{a.tag}] checkpoint unusable ({e}); starting fresh", flush=True)
if resume == 0:
    sim.context.setPositions(pdb.positions)
    print(f"[{a.tag}] {system.getNumParticles()} atoms; minimising", flush=True)
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(310*u.kelvin, a.seed)
    print(f"[{a.tag}] equilibrating 200 ps", flush=True)
    sim.step(int(200/dt_ps))
    sim.currentStep = 0
sim.reporters.append(app.DCDReporter(dcd, int(10/dt_ps), append=os.path.exists(dcd) and resume>0))
sim.reporters.append(app.CheckpointReporter(chk, int(100/dt_ps)))
sim.reporters.append(app.StateDataReporter(
    open(log, "a" if resume>0 else "w"), int(10/dt_ps), step=True, time=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
    speed=True, elapsedTime=True))
remaining = steps_total - resume
print(f"[{a.tag}] running {remaining} steps -> {a.ns} ns total", flush=True)
t0 = time.time()
block = int(1000/dt_ps)   # report every 1 ns
done = resume
while done < steps_total:
    n = min(block, steps_total-done)
    sim.step(n); done += n
    el = time.time()-t0
    print(f"[{a.tag}] {done*dt_ps/1000:7.2f} / {a.ns} ns   "
          f"{(done-resume)*dt_ps/1000/el*86400:6.2f} ns/day   elapsed {el/3600:.2f} h", flush=True)
print(f"[{a.tag}] DONE {done*dt_ps/1000:.2f} ns", flush=True)
