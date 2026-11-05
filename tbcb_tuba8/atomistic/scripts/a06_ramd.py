"""tauRAMD relative residence times (Kokh et al., J Chem Theory Comput 14:3859, 2018) for the alpha tail leaving the
TBCB CAP-Gly groove.

Protocol (GROMACS-RAMD defaults): a force of 585.76 kJ/mol/nm (14 kcal/mol/A) on the tail centre of mass, split
over tail atoms by mass, in a random direction; every 50 steps (100 fs at 2 fs) the tail COM displacement since the
last check is compared with 0.0025 nm and the direction is redrawn if it moved less; a trajectory ends when the
tail COM is 3.0 nm farther from the CAP-Gly COM than at the start, or at 20 ns.  For each system, several
independently equilibrated starts x several RAMD trajectories; tau = median dissociation time over a bootstrap
of trajectories (as in the original tauRAMD analysis), reported with the bootstrap 95 % interval.
These are relative residence times: only ratios between systems are meaningful (koff ratio = 1/tau ratio).

Usage: a06_ramd.py <build_dir> <run_dir> <start> --ntraj 10 --budget 2800 [--tail B --rec A]
"""
import os, sys, json, time, argparse
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

ap = argparse.ArgumentParser()
ap.add_argument('build'); ap.add_argument('run'); ap.add_argument('start', type=int)
ap.add_argument('--ntraj', type=int, default=10); ap.add_argument('--budget', type=float, default=2800.0)
ap.add_argument('--tail', default='B'); ap.add_argument('--rec', default='A')
ap.add_argument('--force', type=float, default=585.76); ap.add_argument('--every', type=int, default=50)
ap.add_argument('--rmin', type=float, default=0.0025); ap.add_argument('--dmax', type=float, default=3.0)
ap.add_argument('--tmax_ns', type=float, default=20.0); ap.add_argument('--platform', default=None)
ap.add_argument('--eq_scale', type=float, default=1.0)
a = ap.parse_args()
t0 = time.time()
W = os.path.join(a.run, f's{a.start}'); os.makedirs(W, exist_ok=True)
log = lambda m: print(m, flush=True)
top, pos = C.load_cif(os.path.join(a.build, 'minimized.cif.gz'))
atoms = list(top.atoms())
lig = [x.index for x in atoms if x.residue.chain.id == a.tail and x.residue.name not in C.SOLVENT]
rec = [x.index for x in atoms if x.residue.chain.id == a.rec and x.residue.name not in C.SOLVENT and x.name == 'CA']
EQ = os.path.join(W, 'equil_state.npz')
seed = 9100 + 37 * a.start
if not os.path.exists(EQ):
    P, V, box, rr = C.equilibrate(top, pos, T=310.0, seed=seed, platform=a.platform, log=log, npt=True, scale=a.eq_scale)
    json.dump(rr, open(os.path.join(W, 'equil.json'), 'w'), indent=1)
    C.save_state_npz(EQ, P, V, box)
P, V, box, _ = C.load_state_npz(EQ)
system = C.make_system(top, hmr=False)
mass = np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in lig])
frac = mass / mass.sum()
f = mm.CustomExternalForce('-frac*(fx*x + fy*y + fz*z)')
for p in ('fx', 'fy', 'fz'):
    f.addGlobalParameter(p, 0.0)
f.addPerParticleParameter('frac')
for i, w in zip(lig, frac):
    f.addParticle(i, [float(w)])
system.addForce(f)
integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
plat, props = C.get_platform(a.platform)
ctx = mm.Context(system, integ, plat, props)
rng = np.random.default_rng(seed)
res_f = os.path.join(W, 'ramd_times.json')
done = json.load(open(res_f)) if os.path.exists(res_f) else []
mrec = np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in rec])
def coms(X):
    return (X[lig] * frac[:, None]).sum(0), (X[rec] * (mrec / mrec.sum())[:, None]).sum(0)
nmax = int(a.tmax_ns * 1000 / 0.002)
while len(done) < a.ntraj and (time.time() - t0) / 60 < a.budget - 10:
    ctx.setPeriodicBoxVectors(*box); ctx.setPositions(P)
    ctx.setVelocitiesToTemperature(310 * unit.kelvin, int(rng.integers(1, 2 ** 30)))
    X = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(C.NM)
    cl, cr = coms(X); d0 = np.linalg.norm(cl - cr); last = cl.copy()
    def newdir():
        v = rng.normal(size=3); v /= np.linalg.norm(v)
        for k, p in zip(v, ('fx', 'fy', 'fz')):
            ctx.setParameter(p, float(k * a.force))
    newdir(); steps = 0; ndir = 1; out = None
    while steps < nmax:
        integ.step(a.every); steps += a.every
        X = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(C.NM)
        cl, cr = coms(X)
        if np.linalg.norm(cl - cr) - d0 > a.dmax:
            out = steps * 0.002 / 1000.0; break
        if np.linalg.norm(cl - last) < a.rmin:
            newdir(); ndir += 1
        last = cl.copy()
    done.append({'t_ns': out, 'censored': out is None, 'n_directions': ndir, 'steps': steps})
    json.dump(done, open(res_f + '.tmp', 'w'), indent=1); os.replace(res_f + '.tmp', res_f)
    log(f'traj {len(done)}: {"%.3f ns" % out if out is not None else "no exit in %.1f ns" % a.tmax_ns}, {ndir} directions')
log(f'{len(done)} trajectories done')
