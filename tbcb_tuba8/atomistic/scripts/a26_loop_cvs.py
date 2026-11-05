"""Loop 32-45 collective variables recomputed from saved trajectories, to calibrate the REST2 state definitions.

Definitions follow a04_rest2.py exactly (native heavy-atom pairs of residues 32-45 with the rest of the alpha chain,
|i-j| > 3, d0 < 4.5 A in the minimised structure, i.e. the coordinates stored in solute.pdb; Best-Hummer-Eaton
switching 1/(1+exp(beta (d - lambda d0))), beta = 50/nm; backbone RMSD of 32-45 after superposing the core CA, i.e.
alpha CA outside 25-55; phi/psi of residue 43), with Q evaluated at lambda = 1.2 (the value a04 uses in-run) and at
lambda = 1.5 and 1.8 (Best, Hummer & Eaton 2013 recommend 1.8 for all-atom models).
Input: a directory holding solute.pdb and trajectories matching --glob (MD: seg_*.dcd; REST2 state 0:
rest2_s0_from*.dcd).  Only the alpha chain (--chain, the letter in solute.pdb) is used.
Output: per-block means (4 contiguous blocks) and the full time series at the chosen stride.
Usage: a26_loop_cvs.py <dir> <out.json> [--glob 'seg_*.dcd'] [--chain A] [--stride 2] [--dt_ps 20]
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md

ap = argparse.ArgumentParser()
ap.add_argument('dir'); ap.add_argument('out')
ap.add_argument('--glob', default='seg_*.dcd'); ap.add_argument('--chain', default='A')
ap.add_argument('--stride', type=int, default=2); ap.add_argument('--dt_ps', type=float, default=20.0)
a = ap.parse_args()
pdb = os.path.join(a.dir, 'solute.pdb')
ch = np.array([l[21] for l in open(pdb) if l.startswith(('ATOM', 'HETATM'))])
ref = md.load(pdb)
top = ref.topology
assert len(ch) == top.n_atoms
fs = sorted(glob.glob(os.path.join(a.dir, a.glob)))
parts = []
for f in fs:
    try:
        parts.append(md.load_dcd(f, top=pdb, stride=a.stride))
    except Exception as ex:                       # a trajectory being copied right now: skip it, say so
        print('skip', f, ex, flush=True)
t = md.join(parts, check_topology=False)
alpha = [x for x in top.atoms if x.residue.is_protein and ch[x.index] == a.chain]
rid = {x.index: x.residue.resSeq for x in alpha}
loop = [x.index for x in alpha if 32 <= x.residue.resSeq <= 45 and x.element.symbol != 'H']
rest = [x.index for x in alpha if not (32 <= x.residue.resSeq <= 45) and x.element.symbol != 'H']
P0 = ref.xyz[0]
d = np.linalg.norm(P0[loop][:, None] - P0[rest][None], axis=2)
ii, jj = np.where(d < 0.45)
keep = [k for k in range(len(ii)) if abs(rid[loop[ii[k]]] - rid[rest[jj[k]]]) > 3]
pairs = np.array([[loop[ii[k]], rest[jj[k]]] for k in keep]); d0 = d[ii[keep], jj[keep]]
D = md.compute_distances(t, pairs, periodic=False)
Q = {lam: (1.0 / (1.0 + np.exp(50.0 * (D - lam * d0)))).mean(axis=1) for lam in (1.2, 1.5, 1.8)}
Q0 = {lam: float(np.mean(1.0 / (1.0 + np.exp(50.0 * (d0 - lam * d0))))) for lam in (1.2, 1.5, 1.8)}
core = [x.index for x in alpha if x.name == 'CA' and not (25 <= x.residue.resSeq <= 55)]
lbb = [x.index for x in alpha if x.name in ('N', 'CA', 'C', 'O') and 32 <= x.residue.resSeq <= 45]
ts = md.Trajectory(t.xyz.copy(), top); ts.superpose(ref, 0, atom_indices=core)
rmsd = np.sqrt(np.mean(np.sum((ts.xyz[:, lbb] - P0[lbb]) ** 2, axis=2), axis=1)) * 10
core_rmsd = md.rmsd(t, ref, 0, atom_indices=core) * 10
by = {(x.residue.resSeq, x.name): x.index for x in alpha}
phi = np.degrees(md.compute_dihedrals(t, np.array([[by[(42, 'C')], by[(43, 'N')], by[(43, 'CA')], by[(43, 'C')]]]))[:, 0])
psi = np.degrees(md.compute_dihedrals(t, np.array([[by[(43, 'N')], by[(43, 'CA')], by[(43, 'C')], by[(44, 'N')]]]))[:, 0])
bl = lambda x: [float(np.mean(b)) for b in np.array_split(np.asarray(x, float), 4)]
r43 = [x for x in alpha if x.residue.resSeq == 43][0].residue.name
out = dict(dir=a.dir, glob=a.glob, files=[os.path.basename(f) for f in fs], n_frames=int(t.n_frames),
           ns_per_frame=a.stride * a.dt_ps / 1000.0, n_native_pairs=int(len(pairs)), res43=r43,
           Q_at_reference=Q0, blocks={'Q_lam1.2': bl(Q[1.2]), 'Q_lam1.5': bl(Q[1.5]), 'Q_lam1.8': bl(Q[1.8]),
                                      'rmsd_loop_A': bl(rmsd), 'rmsd_core_ca_A': bl(core_rmsd),
                                      'p_phi43_pos': bl(phi > 0)},
           means={'Q_lam1.2': float(Q[1.2].mean()), 'Q_lam1.5': float(Q[1.5].mean()), 'Q_lam1.8': float(Q[1.8].mean()),
                  'rmsd_loop_A': float(rmsd.mean()), 'rmsd_core_ca_A': float(core_rmsd.mean()), 'p_phi43_pos': float(np.mean(phi > 0))},
           first_ns={'Q_lam1.2': float(Q[1.2][:max(1, int(1.0 / (a.stride * a.dt_ps / 1000.0)))].mean()),
                     'rmsd_loop_A': float(rmsd[:max(1, int(1.0 / (a.stride * a.dt_ps / 1000.0)))].mean())},
           series={'Q_lam1.2': np.round(Q[1.2], 3).tolist(), 'Q_lam1.8': np.round(Q[1.8], 3).tolist(),
                   'rmsd_loop_A': np.round(rmsd, 2).tolist(), 'phi43': np.round(phi, 1).tolist(), 'psi43': np.round(psi, 1).tolist()})
json.dump(out, open(a.out, 'w'))
print(json.dumps({k: out[k] for k in ('n_frames', 'res43', 'n_native_pairs', 'Q_at_reference', 'means', 'first_ns')}))
