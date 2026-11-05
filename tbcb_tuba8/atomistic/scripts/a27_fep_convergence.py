"""Time-resolved convergence diagnostics for the HREX alchemical legs (reads only fep_ck.npz; no new sampling).

For every replicate of every leg matched by --legs (glob under <runs_root>/fep):
  * the a09 estimate (10 % discard, subsampling by the statistical inefficiency g of the replica-summed own-state
    reduced potential) and pymbar's automated equilibration detection on the same observable (t_eq, g, N_eff);
  * dG in 4 contiguous production blocks, and forward / reverse cumulative estimates (first / last 20,40,60,80,100 %);
  * the free-energy profile split into the physical steps of the 22-state ladder
      new torsions on (0->2), old charges off (2->6), sterics swap (6->15), new charges on (15->19),
      old torsions off (19->21)
    in each half of the production data, so that a half-split drift can be assigned to one step.
Output: one JSON (default <runs_root>/fep_convergence.json) and a short table on stdout.
Usage: a27_fep_convergence.py <runs_root>/fep [--legs '*114*'] [--out x.json] [--T 310]
"""
import os, sys, json, glob, argparse
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument('root'); ap.add_argument('--legs', default='*')
ap.add_argument('--out', default=None); ap.add_argument('--T', type=float, default=310.0)
ap.add_argument('--discard', type=float, default=0.10)
a = ap.parse_args()
a.out = a.out or os.path.join(os.path.dirname(os.path.abspath(a.root.rstrip('/'))), 'fep_convergence.json')
import pymbar
from pymbar import timeseries
kT = 0.0019872043 * a.T

STEPS22 = [('new torsions on', 0, 2), ('old charges off', 2, 6), ('sterics swap', 6, 15),
           ('new charges on', 15, 19), ('old torsions off', 19, 21)]


def analyse(ck):
    z = np.load(ck, allow_pickle=True)
    U = np.array(z['rec_u']); S = np.array(z['rec_s'])
    n, K, _ = U.shape
    if n < 200:
        return None
    t0 = int(a.discard * n)
    own = np.array([U[t, np.arange(K), S[t]].sum() for t in range(n)])
    g = timeseries.statistical_inefficiency(own[t0:])
    stride = max(1, int(np.ceil(g)))
    teq, geq, neff = timeseries.detect_equilibration(own)

    def mbar(t_lo, t_hi):
        ids = np.arange(t_lo, t_hi, stride)
        if len(ids) < 5:
            return None
        u_kn = np.concatenate([U[t].T for t in ids], axis=1)
        st = np.concatenate([S[t] for t in ids])
        N_k = np.bincount(st, minlength=K)
        order = np.argsort(st, kind='stable')
        mb = pymbar.MBAR(u_kn[:, order], N_k, solver_protocol='robust')
        r = mb.compute_free_energy_differences()
        f = r['Delta_f'][0] * kT
        return dict(dG=float(f[-1]), se=float(r['dDelta_f'][0, -1] * kT), f=f.tolist(), n_iter_used=int(len(ids)))

    full = mbar(t0, n)
    edges = np.linspace(t0, n, 5).astype(int)
    blocks = [mbar(edges[i], edges[i + 1]) for i in range(4)]
    fr = (0.2, 0.4, 0.6, 0.8, 1.0)
    fwd = {f: mbar(t0, t0 + int(f * (n - t0))) for f in fr}
    rev = {f: mbar(n - int(f * (n - t0)), n) for f in fr}
    mid = t0 + (n - t0) // 2
    h1, h2 = mbar(t0, mid), mbar(mid, n)
    steps = None
    if K == 22 and h1 and h2:
        steps = {}
        for name, i, j in STEPS22:
            s1 = h1['f'][j] - h1['f'][i]; s2 = h2['f'][j] - h2['f'][i]; sf = full['f'][j] - full['f'][i]
            steps[name] = dict(full=sf, first_half=s1, second_half=s2, change=s2 - s1)
    teq_pymbar = dict(t_eq=int(teq), g=float(geq), N_eff=float(neff))
    return dict(n_iter=int(n), t0=t0, g=float(g), stride=stride, pymbar_equilibration=teq_pymbar,
                dG=full['dG'], se=full['se'],
                eq_detected=mbar(int(teq), n) if teq < n - 5 * stride else None,
                blocks=[None if b is None else dict(dG=b['dG'], se=b['se']) for b in blocks],
                block_edges=edges.tolist(),
                forward={str(f): None if v is None else dict(dG=v['dG'], se=v['se']) for f, v in fwd.items()},
                reverse={str(f): None if v is None else dict(dG=v['dG'], se=v['se']) for f, v in rev.items()},
                halves=[h1['dG'] if h1 else None, h2['dG'] if h2 else None], steps=steps)


res = {}
for leg in sorted(glob.glob(os.path.join(a.root, a.legs))):
    if not os.path.isdir(leg):
        continue
    reps = {}
    for ck in sorted(glob.glob(os.path.join(leg, 'r*', 'fep_ck.npz'))):
        try:
            r = analyse(ck)
        except Exception as ex:
            r = {'error': repr(ex)}
        if r is not None:
            st = os.path.join(os.path.dirname(ck), 'fep_status.json')
            r['status'] = json.load(open(st)) if os.path.exists(st) else {}
            reps[os.path.basename(os.path.dirname(ck))] = r
    res[os.path.basename(leg)] = reps
json.dump(dict(legs=res, T=a.T, discard=a.discard), open(a.out, 'w'), indent=1, default=float)
print('wrote', a.out)
for leg, reps in res.items():
    for rk, r in reps.items():
        if 'error' in r:
            print(leg, rk, 'ERROR', r['error'][:100]); continue
        fmt = lambda d: 'None' if d is None else f"{d['dG']:+.2f}"
        print(f"{leg:30s} {rk} it {r['n_iter']} g {r['g']:.0f} teq {r['pymbar_equilibration']['t_eq']} "
              f"dG {r['dG']:+.3f}+/-{r['se']:.3f} blocks {[fmt(b) for b in r['blocks']]} "
              f"rev {[fmt(r['reverse'][k]) for k in ('0.2', '0.4', '0.6', '0.8', '1.0')]}")
        if r['steps']:
            print('    ' + '  '.join(f"{k}: {v['first_half']:+.2f}->{v['second_half']:+.2f}" for k, v in r['steps'].items()))
