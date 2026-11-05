"""MBAR analysis of the HREX alchemical legs, with the gates the protocol requires.

For each replicate directory (…/r<k>/fep_ck.npz): discard the first 10 % of iterations, subsample by the statistical
inefficiency of the replica-averaged reduced potential, run MBAR (pymbar 4) on all (replica, iteration) samples,
and report dG (state 0 -> K-1) in kcal/mol, the minimum adjacent-state overlap, the forward/reverse
time-split estimates (first half vs second half of the production data), and replica round trips.
Gates: adjacent overlap >= 0.03; |dG(first half) - dG(second half)| <= max(2 * MBAR SE, replicate spread) is
checked at the leg level.  A leg's final value is the mean over replicates that pass, with the SD across them.

Usage: a09_fep_analysis.py <runs_root>/fep  [--T 310] [--out x.json]
  default output: <campaign>/results/fep_results.json when <runs_root>/../results exists (the analysis jobs that were
  submitted before --out was added rely on this), else <runs_root>/fep_results.json
"""
import os, sys, json, glob, argparse
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument('root'); ap.add_argument('--T', type=float, default=310.0)
ap.add_argument('--discard', type=float, default=0.10)
ap.add_argument('--out', default=None, help='output JSON (default: see Usage)')
a = ap.parse_args()
if a.out is None:
    _runs = os.path.dirname(os.path.abspath(a.root.rstrip('/')))
    _res = os.path.join(os.path.dirname(_runs), 'results')
    a.out = os.path.join(_res, 'fep_results.json') if os.path.isdir(_res) else os.path.join(_runs, 'fep_results.json')
import pymbar
from pymbar import timeseries
kT = 0.0019872043 * a.T  # kcal/mol


def analyse(ck):
    z = np.load(ck, allow_pickle=True)
    U = np.array(z['rec_u'])          # (n_iter, K replicas, K states), reduced units
    S = np.array(z['rec_s'])          # (n_iter, K) state of each replica
    n, K, _ = U.shape
    if n < 50:
        return None
    t0 = int(a.discard * n)
    # u of each replica in its own state, summed over replicas, as the decorrelation observable
    own = np.array([U[t, np.arange(K), S[t]].sum() for t in range(n)])
    g = timeseries.statistical_inefficiency(own[t0:])
    idx = np.arange(t0, n, max(1, int(np.ceil(g))))
    def mbar_on(ids):
        u_kn = np.concatenate([U[t].T for t in ids], axis=1)          # (K states, n_samples)
        st = np.concatenate([S[t] for t in ids])
        N_k = np.bincount(st, minlength=K)
        order = np.argsort(st, kind='stable')
        u_kn = u_kn[:, order]
        mb = pymbar.MBAR(u_kn, N_k, solver_protocol='robust')
        r = mb.compute_free_energy_differences()
        ov = mb.compute_overlap()['matrix']
        return float(r['Delta_f'][0, -1]), float(r['dDelta_f'][0, -1]), ov, N_k
    dG, dGse, ov, N_k = mbar_on(idx)
    half = len(idx) // 2
    dG1, _, _, _ = mbar_on(idx[:half]); dG2, _, _, _ = mbar_on(idx[half:])
    adj = [float(ov[k, k + 1]) for k in range(K - 1)]
    # round trips of replicas between state 0 and K-1
    trips = 0
    for r in range(K):
        last = None
        for t in range(n):
            s = S[t, r]
            if s == 0:
                if last == 'top': trips += 1
                last = 'bottom'
            elif s == K - 1 and last == 'bottom':
                last = 'top'
    return dict(n_iter=n, g=float(g), n_used=len(idx), dG_kcal=dG * kT, se_kcal=dGse * kT,
                dG_first_half=dG1 * kT, dG_second_half=dG2 * kT, min_adjacent_overlap=min(adj),
                adjacent_overlap=adj, round_trips=trips, N_k=N_k.tolist())


res = {}
for leg in sorted(glob.glob(os.path.join(a.root, '*'))):
    if not os.path.isdir(leg):
        continue
    reps = {}
    for ck in sorted(glob.glob(os.path.join(leg, 'r*', 'fep_ck.npz'))):
        try:
            r = analyse(ck)
        except Exception as ex:
            r = {'error': str(ex)}
        if r is not None:
            meta = os.path.join(os.path.dirname(ck), 'fep_meta.json')
            st = os.path.join(os.path.dirname(ck), 'fep_status.json')
            r['meta'] = {k: v for k, v in json.load(open(meta)).items() if k != 'states'} if os.path.exists(meta) else {}
            r['status'] = json.load(open(st)) if os.path.exists(st) else {}
            r['gate_overlap'] = r.get('min_adjacent_overlap', 0) >= 0.03
            r['gate_complete'] = bool(r['status'].get('done', False))
            reps[os.path.basename(os.path.dirname(ck))] = r
    ok = [r for r in reps.values() if r.get('gate_overlap') and r.get('gate_complete')]
    vals = np.array([r['dG_kcal'] for r in ok])
    summ = dict(n_replicates=len(reps), n_passing=len(ok),
                mean_dG_kcal=float(vals.mean()) if len(vals) else None,
                sd_dG_kcal=float(vals.std(ddof=1)) if len(vals) > 1 else None,
                sem_dG_kcal=float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else None)
    if len(ok):
        hyst = max(abs(r['dG_first_half'] - r['dG_second_half']) for r in ok)
        summ['max_half_split_kcal'] = hyst
        summ['gate_half_split'] = bool(hyst <= max(2 * max(r['se_kcal'] for r in ok), summ['sd_dG_kcal'] or 0))
    res[os.path.basename(leg)] = dict(summary=summ, replicates=reps)
# ddG_fold = dG(folded) - dG(unfolded) for every mutation with both legs
dd = {}
for k, v in res.items():
    if k.endswith('_folded'):
        base = k[:-len('_folded')]
        mut = base.split('_')[-1]
        unf = [u for u in res if u.endswith(f'{mut}_unfolded')]
        f = v['summary']
        for u in unf:
            g = res[u]['summary']
            if f['mean_dG_kcal'] is not None and g['mean_dG_kcal'] is not None:
                sd = np.sqrt((f['sd_dG_kcal'] or 0) ** 2 + (g['sd_dG_kcal'] or 0) ** 2)
                dd[f'{base}  [ref {u}]'] = dict(ddG_fold_kcal=f['mean_dG_kcal'] - g['mean_dG_kcal'], sd_kcal=float(sd),
                                               folded=f, unfolded=g)
out = dict(legs=res, ddG_fold=dd, T=a.T, note='ddG_fold > 0 = mutation destabilises the folded state')
json.dump(out, open(a.out, 'w'), indent=1, default=float)
print('wrote', a.out)
for k, v in dd.items():
    print(f"{k}: ddG_fold = {v['ddG_fold_kcal']:+.2f} +/- {v['sd_kcal']:.2f} kcal/mol")
for k, v in res.items():
    s = v['summary']
    print(k, s)
