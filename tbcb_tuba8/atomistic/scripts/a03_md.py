"""Unbiased explicit-solvent MD of one build: staged equilibration, then production in resumable segments.

Usage: a03_md.py <build_dir> <run_dir> <replica> <target_ns> <budget_min> [--T 310] [--nvt] [--save_ps 20]
  310 K NPT production is the unbiased ensemble used for the MSM and for local-unfolding statistics.
  --T 450 --nvt gives the high-temperature unfolding runs (seeds and unfolding order only, not a physical path).
Outputs in <run_dir>/r<replica>/: equil.json, state.npz (latest checkpoint), seg_XXX.dcd (solute + Mg-shell waters),
  seg_XXX_mg.csv (Mg2+ to its starting shell oxygens), progress.json, solute.pdb (topology for the DCDs).
"""
import os, sys, json, time, argparse, glob, shutil
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

ap = argparse.ArgumentParser()
ap.add_argument('build'); ap.add_argument('run'); ap.add_argument('rep', type=int)
ap.add_argument('target_ns', type=float); ap.add_argument('budget_min', type=float, nargs='?', default=None)
ap.add_argument('--budget', type=float, default=None, help='wall-clock budget in minutes (same as the positional)')
ap.add_argument('--T', type=float, default=310.0); ap.add_argument('--nvt', action='store_true')
ap.add_argument('--save_ps', type=float, default=20.0); ap.add_argument('--platform', default=None)
ap.add_argument('--start', default=None, help='optional start structure (cif.gz) instead of the build')
ap.add_argument('--eq_scale', type=float, default=1.0)
a = ap.parse_args()
a.budget_min = a.budget if a.budget is not None else (a.budget_min if a.budget_min is not None else 2800.0)
t0 = time.time()
W = os.path.join(a.run, f'r{a.rep}'); os.makedirs(W, exist_ok=True)
LOCAL = os.path.join(os.environ.get('L_SCRATCH_JOB') or os.environ.get('TMPDIR') or '/tmp', f'md_{os.path.basename(a.run)}_r{a.rep}')
os.makedirs(LOCAL, exist_ok=True)
log = lambda m: print(m, flush=True)
start = a.start or (os.path.join(a.build, 'minimized.cif.gz') if os.path.exists(os.path.join(a.build, 'minimized.cif.gz'))
                    else os.path.join(a.build, 'solvated.cif.gz'))
top, pos = C.load_cif(start)
atoms = list(top.atoms())
shell = C.mg_shell_pairs(top, pos)
save_idx = sorted(set(C.solute_indices(top)) | {b for _, b in shell} |
                  {x.index for _, b in shell for x in atoms[b].residue.atoms()})
if not os.path.exists(os.path.join(W, 'solute.pdb')):
    m = app.Modeller(top, pos)
    keep = set(save_idx)
    m.delete([x for x in m.topology.atoms() if x.index not in keep])
    with open(os.path.join(W, 'solute.pdb'), 'w') as fh:
        app.PDBFile.writeFile(m.topology, m.positions, fh)
sub_top = app.PDBFile(os.path.join(W, 'solute.pdb')).topology
CK = os.path.join(W, 'state.npz'); LCK = os.path.join(LOCAL, 'state.npz')
prog_f = os.path.join(W, 'progress.json')
prog = json.load(open(prog_f)) if os.path.exists(prog_f) else {'ns_done': 0.0, 'segments': 0}
seed = 1000 + 17 * a.rep + int(a.T)
if not os.path.exists(CK):
    log(f'equilibrating {start} at {a.T if not a.nvt else 310.0} K, seed {seed}')
    P, V, box, rec = C.equilibrate(top, pos, T=310.0, seed=seed, platform=a.platform, log=log, npt=True, scale=a.eq_scale)
    rec['start'] = start
    json.dump(rec, open(os.path.join(W, 'equil.json'), 'w'), indent=1)
    C.save_state_npz(LCK, P, V, box, extra={'ns_done': 0.0})
    C.push(LCK, CK)
P, V, box, extra = C.load_state_npz(CK)
ns_done = float(extra.get('ns_done', prog['ns_done']))
if ns_done >= a.target_ns:
    log(f'already at {ns_done} ns >= target'); sys.exit(0)
system = C.make_system(top, hmr=True)
if not a.nvt:
    system.addForce(mm.MonteCarloBarostat(1 * unit.bar, a.T * unit.kelvin, 25))
dt = 0.004 if a.T <= 330 else 0.002
integ = mm.LangevinMiddleIntegrator(a.T * unit.kelvin, 1 / unit.picosecond, dt * unit.picoseconds)
integ.setRandomNumberSeed(seed + prog['segments'] * 7919)
plat, props = C.get_platform(a.platform)
ctx = mm.Context(system, integ, plat, props)
ctx.setPeriodicBoxVectors(*box); ctx.setPositions(P)
if prog['segments'] == 0 and a.T != 310.0:
    ctx.setVelocitiesToTemperature(a.T * unit.kelvin, seed)
else:
    ctx.setVelocities(V)
dev = plat.getPropertyValue(ctx, 'DeviceName') if plat.getName() in ('CUDA', 'OpenCL') else 'cpu'
seg = prog['segments']
# claim the segment number before writing anything, so a job killed mid-segment can never be overwritten
prog['segments'] = seg + 1
prog.setdefault('seg_start_ns', {})[str(seg)] = ns_done
json.dump(prog, open(prog_f + '.tmp', 'w'), indent=1); os.replace(prog_f + '.tmp', prog_f)
dcd_local = os.path.join(LOCAL, f'seg_{seg:03d}.dcd'); mg_local = os.path.join(LOCAL, f'seg_{seg:03d}_mg.csv')
dcd = app.DCDFile(open(dcd_local, 'wb'), sub_top, dt * int(round(a.save_ps / dt)))
mgf = open(mg_local, 'w'); mgf.write('t_ns,' + ','.join(f'{atoms[j].residue.name}{atoms[j].residue.id}:{atoms[j].name}' for _, j in shell) + ',n_O_within_2.5A\n')
nsave = int(round(a.save_ps / dt))
chunk = nsave * max(1, int(round(100.0 / a.save_ps)))
allO = np.array([x.index for x in atoms if x.element is not None and x.element.symbol == 'O'])
mg_idx = [i for i, _ in shell][:1]
log(f'production on {plat.getName()} {dev}: {ns_done:.2f} -> {a.target_ns} ns, dt {dt} ps, T {a.T}, {"NVT" if a.nvt else "NPT"}')
tlast = time.time(); tw = time.time(); steps_done_here = 0
while ns_done < a.target_ns - 1e-9:
    for _ in range(chunk // nsave):
        integ.step(nsave); steps_done_here += nsave
        ns_done += nsave * dt / 1000.0
        st = ctx.getState(getPositions=True, enforcePeriodicBox=False)
        X = st.getPositions(asNumpy=True).value_in_unit(C.NM)
        dcd.writeModel(X[save_idx] * C.NM, periodicBoxVectors=st.getPeriodicBoxVectors())
        if mg_idx:
            m0 = X[mg_idx[0]]
            bx = np.array([v.value_in_unit(C.NM) for v in st.getPeriodicBoxVectors()])
            d = [np.linalg.norm(X[j] - m0) * 10 for _, j in shell]
            dO = X[allO] - m0
            # minimum image for a triclinic box (reduced form, as OpenMM uses)
            for k in (2, 1, 0):
                dO -= np.outer(np.round(dO[:, k] / bx[k, k]), bx[k])
            nO = int(np.sum(np.linalg.norm(dO, axis=1) < 0.25))
            mgf.write(f'{ns_done:.4f},' + ','.join(f'{x:.3f}' for x in d) + f',{nO}\n')
        if ns_done >= a.target_ns - 1e-9:
            break
    if not np.all(np.isfinite(X)):
        raise RuntimeError('non-finite coordinates')
    if time.time() - tlast > 600 or ns_done >= a.target_ns - 1e-9 or (time.time() - t0) / 60 > a.budget_min - 4:
        st = ctx.getState(getPositions=True, getVelocities=True)
        C.save_state_npz(LCK, st.getPositions(asNumpy=True), st.getVelocities(asNumpy=True), st.getPeriodicBoxVectors(),
                         extra={'ns_done': ns_done})
        mgf.flush(); dcd._file.flush() if hasattr(dcd, '_file') else None
        C.push(LCK, CK); C.push(dcd_local, os.path.join(W, os.path.basename(dcd_local)))
        C.push(mg_local, os.path.join(W, os.path.basename(mg_local)))
        rate = steps_done_here * dt / 1000.0 / ((time.time() - tw) / 86400)
        prog.update(ns_done=ns_done, ns_per_day=round(rate, 1), device=dev, updated=time.ctime())
        prog.setdefault('seg_end_ns', {})[str(seg)] = ns_done
        json.dump(prog, open(prog_f + '.tmp', 'w'), indent=1); os.replace(prog_f + '.tmp', prog_f)
        log(f'{ns_done:.2f} ns  {rate:.0f} ns/day')
        tlast = time.time()
        if (time.time() - t0) / 60 > a.budget_min - 4:
            break
mgf.close()
json.dump(prog, open(prog_f + '.tmp', 'w'), indent=1); os.replace(prog_f + '.tmp', prog_f)
C.push(dcd_local, os.path.join(W, os.path.basename(dcd_local))); C.push(mg_local, os.path.join(W, os.path.basename(mg_local)))
log(f'segment {seg} done at {ns_done:.2f} ns')
