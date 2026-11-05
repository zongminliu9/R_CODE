"""GPU sanity check: platform, openmmtools import, ns/day for one build with HMR 4 fs."""
import os, sys, time, gzip
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a02_minimize import forcefield, load
print('openmm', mm.__version__, [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())], flush=True)
try:
    import openmmtools, pymbar
    print('openmmtools', openmmtools.__version__, 'pymbar', pymbar.__version__)
except Exception as ex:
    print('openmmtools import failed', ex)
top, pos = load(sys.argv[1])
ff = forcefield()
s = ff.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds,
                    hydrogenMass=3.0 * unit.amu)
s.addForce(mm.MonteCarloBarostat(1 * unit.bar, 310 * unit.kelvin, 25))
integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds)
plat = mm.Platform.getPlatformByName(sys.argv[2] if len(sys.argv) > 2 else 'CUDA')
sim = app.Simulation(top, s, integ, plat, {'Precision': 'mixed'} if plat.getName() in ('CUDA', 'OpenCL') else {})
sim.context.setPositions(pos); sim.context.setPeriodicBoxVectors(*top.getPeriodicBoxVectors())
mm.LocalEnergyMinimizer.minimize(sim.context, 10.0, 200)
sim.context.setVelocitiesToTemperature(310 * unit.kelvin)
sim.step(500)
t = time.time(); n = 5000; sim.step(n); dt = time.time() - t
print(f'{s.getNumParticles()} atoms  {n * 0.004 / 1000 / dt * 86400:.1f} ns/day on {plat.getName()} '
      f'({plat.getPropertyValue(sim.context, "DeviceName") if plat.getName() in ("CUDA", "OpenCL") else ""})', flush=True)
