"""High-temperature (450 K, NVT) unfolding trajectories: which parts of alpha-tubulin lose structure first.
These runs are not a physical folding pathway; they rank how easily each region gives way and provide seeds.
Per trajectory and per folding unit (R_CODE results/structure/folding_units.json, 6 spectral units; residues > 437
dropped): fraction of native heavy-atom contacts Q_u(t) (contacts in the minimised start, |i-j| > 3, 4.5 A, Best-
Hummer switching), helix/strand content (DSSP), and t_half(u) = first time the 1-ns running mean of Q_u falls below 0.5.
Also the 32-45 loop, the 110-120 segment (I114 environment), and GTP contacts.  Unfolding order = units sorted by
t_half; a unit that never falls below 0.5 within the run is censored at the run length.
Usage: a18_unfold_analysis.py <md450_root>/<sid> <build_dir> <out.json>
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md

ap = argparse.ArgumentParser(); ap.add_argument('mdsys'); ap.add_argument('build'); ap.add_argument('out')
ap.add_argument('--dt_ps', type=float, default=10.0)
a = ap.parse_args()
FU = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'inputs', 'folding_units.json')))
unit_of = {r: u for r, u in zip(FU['residues'], FU['unit']) if r <= 437}
groups = {f'unit{u}': sorted(r for r, x in unit_of.items() if x == u) for u in sorted(set(unit_of.values()))}
groups['loop32_45'] = list(range(32, 46)); groups['seg110_120'] = list(range(110, 121))
out = {'groups': {k: [min(v), max(v), len(v)] for k, v in groups.items()}, 'replicas': {}}
for rdir in sorted(glob.glob(os.path.join(a.mdsys, 'r*'))):
    segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd'))); pdb = os.path.join(rdir, 'solute.pdb')
    if not segs:
        continue
    t = md.join([md.load_dcd(s, top=pdb) for s in segs], check_topology=False)
    ref = md.load(pdb)
    top = t.topology
    heavy = top.select('protein and not element H')
    rs = np.array([top.atom(i).residue.resSeq for i in heavy])
    d0 = md.compute_distances(ref, np.array([[i, j] for i in heavy[::1] for j in []]).reshape(0, 2)) if False else None
    # native contacts per group: pairs (i in group, j anywhere) within 4.5 A in the reference
    X0 = ref.xyz[0][heavy]
    rec = {'n_frames': int(t.n_frames), 'ns': float(t.n_frames * a.dt_ps / 1000)}
    for g, res in groups.items():
        gi = np.where(np.isin(rs, res))[0]
        if len(gi) == 0:
            continue
        dd = np.linalg.norm(X0[gi][:, None] - X0[None], axis=2)
        ii, jj = np.where((dd < 0.45) & (np.abs(rs[gi][:, None] - rs[None]) > 3))
        pairs = np.stack([heavy[gi][ii], heavy[jj]], axis=1); r0 = dd[ii, jj]
        if len(pairs) == 0:
            continue
        D = md.compute_distances(t, pairs, periodic=False)
        Q = np.mean(1.0 / (1.0 + np.exp(50.0 * (D - 1.2 * r0[None]))), axis=1)
        w = max(1, int(1000 / a.dt_ps))
        rm = np.convolve(Q, np.ones(w) / w, mode='valid')
        below = np.where(rm < 0.5)[0]
        thalf = float((below[0] + w / 2) * a.dt_ps / 1000) if len(below) else None
        rec[g] = {'n_native_contacts': int(len(pairs)), 'Q_final_1ns': float(Q[-w:].mean()), 't_half_ns': thalf}
    ss = md.compute_dssp(t, simplified=True)
    resseq = np.array([r.resSeq for r in top.residues if r.is_protein])
    for g, res in groups.items():
        m = np.isin(resseq, res)
        if m.any() and g in rec:
            rec[g]['helix_start'] = float(np.mean(ss[:10, m] == 'H')); rec[g]['helix_end'] = float(np.mean(ss[-max(1, int(1000 / a.dt_ps)):, m] == 'H'))
            rec[g]['strand_start'] = float(np.mean(ss[:10, m] == 'E')); rec[g]['strand_end'] = float(np.mean(ss[-max(1, int(1000 / a.dt_ps)):, m] == 'E'))
    gtp = top.select('resname GTP and not element H'); prot = top.select('protein and not element H')
    if len(gtp):
        D = md.compute_neighbors(t[::10], 0.45, gtp, haystack_indices=prot)
        rec['gtp_protein_contacts_first_last'] = [len(D[0]), len(D[-1])]
    order = sorted([g for g in groups if g in rec and g.startswith('unit')],
                   key=lambda g: rec[g]['t_half_ns'] if rec[g]['t_half_ns'] is not None else 1e9)
    rec['unfolding_order_units'] = order
    out['replicas'][os.path.basename(rdir)] = rec
    print(os.path.basename(rdir), order, {g: rec[g]['t_half_ns'] for g in order}, flush=True)
json.dump(out, open(a.out, 'w'), indent=1)
