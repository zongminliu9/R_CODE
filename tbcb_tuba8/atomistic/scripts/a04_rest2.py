"""REST2 (Wang, Friesner & Berne, J Phys Chem B 115:9431, 2011) on a 'hot' region of the alpha-tubulin monomer, with
GTP.Mg2+ bound, explicit TIP3P water kept at 310 K.

Scaling, for hot atoms only: charges x sqrt(lam), LJ epsilon x lam, so hot-hot nonbonded scale by lam and hot-cold
by sqrt(lam); 1-4 exceptions inside the region x lam, across its boundary x sqrt(lam); proper/improper torsions x lam
when both central atoms are hot and x sqrt(lam) when one is.  Bonds and angles are not scaled.  Effective solute
temperature T_eff = 310 K / lam.  State 0 (lam = 1) is the physical Hamiltonian.  NVT at the equilibrated box;
the dispersion correction is switched off so that it cannot depend on lam.

Usage: a04_rest2.py <build_dir> <run_dir> <rep> --hot 30-50 --n 12 --tmax 500 --ns 30 --budget 2800
Records per iteration (2 ps): the reduced potential of every replica in every state (MBAR), the state permutation,
and per-replica CVs: Q of the 32-45 loop (native contacts with the rest of the protein), Q inside the hot region,
backbone RMSD of 32-45 after superposition on the core, and phi/psi of residues 34, 40, 43, 44, 45.
"""
import os, sys, json, time, argparse, math
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from hrex import HREX

ap = argparse.ArgumentParser()
ap.add_argument('build'); ap.add_argument('run'); ap.add_argument('rep', type=int)
ap.add_argument('--hot', default='30-50'); ap.add_argument('--n', type=int, default=12)
ap.add_argument('--tmax', type=float, default=500.0); ap.add_argument('--ns', type=float, default=30.0)
ap.add_argument('--budget', type=float, default=2800.0); ap.add_argument('--steps', type=int, default=500)
ap.add_argument('--platform', default=None); ap.add_argument('--eq_scale', type=float, default=1.0)
a = ap.parse_args()
W = os.path.join(a.run, f'r{a.rep}'); os.makedirs(W, exist_ok=True)
LOCAL = os.path.join(os.environ.get('L_SCRATCH_JOB') or os.environ.get('TMPDIR') or '/tmp', f'rest2_{os.path.basename(a.run)}_r{a.rep}')
os.makedirs(LOCAL, exist_ok=True)
log = lambda m: print(m, flush=True)
lo, hi = [int(x) for x in a.hot.split('-')]
start = os.path.join(a.build, 'minimized.cif.gz')
top, pos = C.load_cif(start)
atoms = list(top.atoms())
prot_res = [r for r in top.residues() if r.name not in C.SOLVENT | {'GTP', 'MG'}]
resid = lambda r: int(r.id)
hot = [x.index for r in prot_res if lo <= resid(r) <= hi for x in r.atoms()]
hotset = set(hot)
EQ = os.path.join(W, 'equil_state.npz')
seed = 4242 + 101 * a.rep
if not os.path.exists(EQ):
    P, V, box, rec = C.equilibrate(top, pos, T=310.0, seed=seed, platform=a.platform, log=log, npt=True,
                                   scale=a.eq_scale)
    json.dump(rec, open(os.path.join(W, 'equil.json'), 'w'), indent=1)
    C.save_state_npz(EQ, P, V, box)
P, V, box, _ = C.load_state_npz(EQ)

system = C.make_system(top, hmr=True)
nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
nb.setUseDispersionCorrection(False)
nb.addGlobalParameter('rest_sq', 1.0); nb.addGlobalParameter('rest_lam', 1.0)
for i in hot:
    q, s, e = nb.getParticleParameters(i)
    nb.setParticleParameters(i, 0.0, s, 0.0)
    nb.addParticleParameterOffset('rest_sq', i, q, 0.0, 0.0)
    nb.addParticleParameterOffset('rest_lam', i, 0.0, 0.0, e)
nexc = [0, 0]
for k in range(nb.getNumExceptions()):
    i, j, qq, s, e = nb.getExceptionParameters(k)
    hi_, hj = i in hotset, j in hotset
    if not (hi_ or hj):
        continue
    if qq.value_in_unit(unit.elementary_charge ** 2) == 0 and e.value_in_unit(unit.kilojoule_per_mole) == 0:
        continue
    par = 'rest_lam' if (hi_ and hj) else 'rest_sq'
    nb.setExceptionParameters(k, i, j, 0.0, s, 0.0)
    nb.addExceptionParameterOffset(par, k, qq, 0.0, e)
    nexc[0 if par == 'rest_lam' else 1] += 1
pt = [f for f in system.getForces() if isinstance(f, mm.PeriodicTorsionForce)][0]
ct = mm.CustomTorsionForce('(sel*rest_sq + (1-sel)*rest_lam)*kt*(1+cos(per*theta-ph))')
ct.addGlobalParameter('rest_sq', 1.0); ct.addGlobalParameter('rest_lam', 1.0)
for p in ('sel', 'kt', 'per', 'ph'):
    ct.addPerTorsionParameter(p)
ntor = [0, 0]
for k in range(pt.getNumTorsions()):
    i, j, l, m, per, ph, kk = pt.getTorsionParameters(k)
    c = (j in hotset) + (l in hotset)
    if c == 0:
        continue
    ct.addTorsion(i, j, l, m, [1.0 if c == 1 else 0.0, kk.value_in_unit(unit.kilojoule_per_mole), float(per),
                               ph.value_in_unit(unit.radian)])
    pt.setTorsionParameters(k, i, j, l, m, per, ph, 0.0)
    ntor[c - 1] += 1
system.addForce(ct)
Teff = np.geomspace(310.0, a.tmax, a.n)
states = [{'rest_lam': 310.0 / T, 'rest_sq': math.sqrt(310.0 / T)} for T in Teff]

# ---------------- CVs
prot_atoms = [x for r in prot_res for x in r.atoms()]
heavy = [x.index for x in prot_atoms if x.element.symbol != 'H']
P0 = np.array(pos.value_in_unit(C.NM))
loop = [x.index for x in prot_atoms if 32 <= int(x.residue.id) <= 45 and x.element.symbol != 'H']
rest_h = [i for i in heavy if not (32 <= int(atoms[i].residue.id) <= 45)]
hot_h = [i for i in hot if atoms[i].element.symbol != 'H']
def native_pairs(A, B, same=False):
    XA, XB = P0[A], P0[B]
    d = np.linalg.norm(XA[:, None] - XB[None], axis=2)
    out = []
    for ii, jj in zip(*np.where(d < 0.45)):
        ra, rb = int(atoms[A[ii]].residue.id), int(atoms[B[jj]].residue.id)
        if abs(ra - rb) > 3 and (not same or A[ii] < B[jj]):
            out.append((A[ii], B[jj], d[ii, jj]))
    return out
nat_loop = native_pairs(loop, rest_h)
nat_hot = native_pairs(hot_h, hot_h, same=True)
core_ca = [x.index for x in prot_atoms if x.name == 'CA' and not (25 <= int(x.residue.id) <= 55)]
loop_bb = [x.index for x in prot_atoms if x.name in ('N', 'CA', 'C', 'O') and 32 <= int(x.residue.id) <= 45]
byres = {(int(x.residue.id), x.name): x.index for x in prot_atoms}
TORS = []
for r in (34, 40, 43, 44, 45):
    TORS.append((f'phi{r}', [byres[(r - 1, 'C')], byres[(r, 'N')], byres[(r, 'CA')], byres[(r, 'C')]]))
    TORS.append((f'psi{r}', [byres[(r, 'N')], byres[(r, 'CA')], byres[(r, 'C')], byres[(r + 1, 'N')]]))
def Q(X, pairs):
    if not pairs:
        return float('nan')
    i = np.array([p[0] for p in pairs]); j = np.array([p[1] for p in pairs]); d0 = np.array([p[2] for p in pairs])
    d = np.linalg.norm(X[i] - X[j], axis=1)
    return float(np.mean(1.0 / (1.0 + np.exp(50.0 * (d - 1.2 * d0)))))   # Best-Hummer-Eaton switching, beta = 50/nm
ref_core = P0[core_ca]; ref_loop = P0[loop_bb]
def cvs(X):
    A = X[core_ca]; ca, cb = A.mean(0), ref_core.mean(0)
    U, S, Vt = np.linalg.svd((A - ca).T @ (ref_core - cb)); dd = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, dd]) @ U.T
    L = (R @ (X[loop_bb] - ca).T).T + cb
    out = [Q(X, nat_loop), Q(X, nat_hot), float(np.sqrt(np.mean(np.sum((L - ref_loop) ** 2, 1)))) * 10]
    out += [C.dihedral(*[X[k] for k in idx]) for _, idx in TORS]
    return out
CV_NAMES = ['Q_loop', 'Q_hot', 'rmsd_loop_A'] + [n for n, _ in TORS]
meta = dict(build=a.build, hot=a.hot, n_hot_atoms=len(hot), Teff=Teff.tolist(), lam=[s['rest_lam'] for s in states],
            exceptions_scaled_lam_sqrt=nexc, torsions_scaled_lam_sqrt=ntor, n_native_loop=len(nat_loop),
            n_native_hot=len(nat_hot), cv_names=CV_NAMES, steps_per_iter=a.steps, dt_ps=0.004)
json.dump(meta, open(os.path.join(W, 'rest2_meta.json'), 'w'), indent=1)
log(json.dumps(meta))
solute = [x.index for x in prot_atoms] + [x.index for x in atoms if x.residue.name in ('GTP', 'MG')]
mm_ = app.Modeller(top, pos); keep = set(solute)
mm_.delete([x for x in mm_.topology.atoms() if x.index not in keep])
if not os.path.exists(os.path.join(W, 'solute.pdb')):
    with open(os.path.join(W, 'solute.pdb'), 'w') as fh:
        app.PDBFile.writeFile(mm_.topology, mm_.positions, fh)
n_iter = int(round(a.ns * 1000 / (a.steps * 0.004)))
h = HREX(system, P, box, states, temperature_K=310.0, timestep_ps=0.004, steps_per_iter=a.steps, platform=a.platform,
         shared_dir=W, local_dir=LOCAL, tag='rest2', cv_fn=cvs, save_states=(0, a.n - 1), save_idx=np.array(solute),
         save_topology=mm_.topology, seed=seed)
h.resume_or_init()
if h.iteration == 0:
    h.vel[:] = V  # start every replica from the equilibrated velocities; the thermostat decorrelates them
done = h.run(n_iter, a.budget, ck_every_min=10.0, dcd_every=1, log_every=100)
json.dump({'iteration': h.iteration, 'target': n_iter, 'done': bool(done), 'acceptance': h.acceptance(),
           'platform': h.platform, 'time': time.ctime()}, open(os.path.join(W, 'rest2_status.json'), 'w'), indent=1)
log(f'iteration {h.iteration}/{n_iter} done={done} acceptance {np.round(h.acceptance(), 2).tolist()}')
