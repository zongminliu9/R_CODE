"""Representative alpha-tubulin intermediates for ensemble docking.

From REST2 (state-0 = physical 310 K trajectory, rest2_s0_*.dcd, paired with the per-iteration CV record): frames are
classed N / I / U by the loop state definition of a10 (Q_loop, loop RMSD) and, within each class, clustered by loop
backbone RMSD after core superposition (k-medoids on a subsample); the medoid of each of the most populated clusters
is written (protein + GTP + Mg, PDB), with its population share.
From the 450 K runs (optional): the frame where folding unit 2 (contains the 32-45 loop and 17 divergent positions)
first has Q < 0.5 while the rest keeps Q > 0.7 - a partially unfolded N-domain intermediate seed.
Usage: a23_representatives.py <rest2_run/r0> <out_dir> [--per_class 2]
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md

ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('out'); ap.add_argument('--per_class', type=int, default=2)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
z = np.load(os.path.join(a.run, 'rest2_ck.npz'), allow_pickle=True)
meta = json.load(open(os.path.join(a.run, 'rest2_meta.json')))
S = np.array(z['rec_s']); CV = np.array([np.array(x, float) for x in z['rec_cv']])
names = meta['cv_names']; ix = {k: i for i, k in enumerate(names)}
n = len(S)
# CV of the replica at state 0 for every iteration, in order: this matches the state-0 DCD frame order
cv0 = np.array([CV[t][np.where(S[t] == 0)[0][0]] for t in range(n)])
dcds = sorted(glob.glob(os.path.join(a.run, 'rest2_s0_from*.dcd')))
pdb = os.path.join(a.run, 'solute.pdb')
traj = md.join([md.load_dcd(f, top=pdb) for f in dcds], check_topology=False)
m = min(len(traj), len(cv0)); traj = traj[:m]; cv0 = cv0[:m]
start = int(0.1 * m)
q, r = cv0[:, ix['Q_loop']], cv0[:, ix['rmsd_loop_A']]
lab = np.where((q >= 0.75) & (r < 2.5), 'N', np.where(q < 0.45, 'U', 'I'))
top = traj.topology
core = top.select('name CA and (resSeq < 25 or resSeq > 55)'); loop = top.select('backbone and resSeq 32 to 45')
traj.superpose(traj, 0, atom_indices=core)
reps = []
for c in ('N', 'I', 'U'):
    idx = np.where(lab[start:] == c)[0] + start
    if len(idx) < 5:
        continue
    sub = idx[:: max(1, len(idx) // 400)]
    X = traj.xyz[sub][:, loop]
    D = np.sqrt(np.mean(np.sum((X[:, None] - X[None]) ** 2, axis=3), axis=2)) * 10
    # k-medoids (simple PAM-like iterations)
    k = min(a.per_class, len(sub)); rng = np.random.default_rng(0)
    med = list(rng.choice(len(sub), k, replace=False))
    for _ in range(30):
        assign = np.argmin(D[:, med], axis=1)
        new = [int(np.argmin(D[np.ix_(assign == j, assign == j)].sum(axis=1)) if np.any(assign == j) else med[j]) for j in range(k)]
        new = [np.where(assign == j)[0][new[j]] if np.any(assign == j) else med[j] for j in range(k)]
        if new == med:
            break
        med = new
    assign = np.argmin(D[:, med], axis=1)
    for j, mi in enumerate(med):
        fr = int(sub[mi]); share = float(np.mean(assign == j)) * len(idx) / (m - start)
        fn = os.path.join(a.out, f'rep_{c}{j}.pdb')
        traj[fr].save_pdb(fn)
        reps.append(dict(file=os.path.basename(fn), cls=c, frame=fr, share_of_310K_frames=share,
                         Q_loop=float(q[fr]), rmsd_loop_A=float(r[fr]), phi43=float(cv0[fr, ix['phi43']])))
json.dump(dict(run=a.run, n_frames=int(m), class_fraction={c: float(np.mean(lab[start:] == c)) for c in ('N', 'I', 'U')},
               representatives=reps), open(os.path.join(a.out, 'representatives.json'), 'w'), indent=1)
print(json.dumps(reps, indent=1))
