"""Analyse a REST2 run: populations and free energies of H1-S2 loop states at 310 K, from MBAR over every replica.

Samples: every (iteration, replica) after a 10 % burn-in, subsampled by the statistical inefficiency of the summed
reduced potential.  MBAR weights for state 0 (the physical 310 K Hamiltonian) reweight all samples.
Loop states (defined on CVs computed in-run against the minimised starting structure):
  N  native-like   : Q_loop >= 0.75 and loop backbone RMSD < 2.5 A
  U  loop-unfolded : Q_loop < 0.45
  I  intermediate  : everything else
and, separately, residue 43 handedness (phi43 > 0, left-handed / Gly-like, versus phi43 < 0).
dG_s = -kT ln(p_s / p_N).  Uncertainty: SD over 4 contiguous blocks of iterations (MBAR redone per block).
PMFs along Q_loop and phi43 (kcal/mol), and the barrier between the N and U basins along Q_loop when both are seen.
Usage: a10_rest2_analysis.py <rest2 run dir (…/r0)> [--T 310]
"""
import os, sys, json, argparse
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('--T', type=float, default=310.0)
ap.add_argument('--discard', type=float, default=0.1)
a = ap.parse_args()
import pymbar
from pymbar import timeseries
kT = 0.0019872043 * a.T
z = np.load(os.path.join(a.run, 'rest2_ck.npz'), allow_pickle=True)
meta = json.load(open(os.path.join(a.run, 'rest2_meta.json')))
U = np.array(z['rec_u']); S = np.array(z['rec_s']); CV = np.array([np.array(x, float) for x in z['rec_cv']])
n, K, _ = U.shape
names = meta['cv_names']; ix = {k: i for i, k in enumerate(names)}
t0 = int(a.discard * n)
own = np.array([U[t, np.arange(K), S[t]].sum() for t in range(n)])
g = float(timeseries.statistical_inefficiency(own[t0:])) if n - t0 > 20 else 1.0
iters = np.arange(t0, n, max(1, int(np.ceil(g))))


def mbar_w(its):
    u_kn = np.concatenate([U[t].T for t in its], axis=1)
    st = np.concatenate([S[t] for t in its])
    N_k = np.bincount(st, minlength=K)
    mb = pymbar.MBAR(u_kn, N_k, solver_protocol='robust')
    w = mb.weights()[:, 0]
    cv = np.concatenate([CV[t] for t in its], axis=0)
    return w / w.sum(), cv, mb


def states(cv):
    q, r, phi = cv[:, ix['Q_loop']], cv[:, ix['rmsd_loop_A']], cv[:, ix['phi43']]
    lab = np.full(len(q), 'I', dtype=object)
    lab[(q >= 0.75) & (r < 2.5)] = 'N'; lab[q < 0.45] = 'U'
    return lab, phi > 0


def summary(w, cv):
    lab, left = states(cv)
    p = {s: float(w[lab == s].sum()) for s in ('N', 'I', 'U')}
    out = {'p': p, 'p_phi43_left': float(w[left].sum()),
           'mean_Q_loop': float(np.sum(w * cv[:, ix['Q_loop']])), 'mean_rmsd_loop_A': float(np.sum(w * cv[:, ix['rmsd_loop_A']]))}
    out['dG_kcal'] = {s: (-kT * np.log(p[s] / p['N']) if p[s] > 0 and p['N'] > 0 else None) for s in ('I', 'U')}
    pl = out['p_phi43_left']
    out['dG_phi43_left_vs_right_kcal'] = -kT * np.log(pl / (1 - pl)) if 0 < pl < 1 else None
    return out


w, cv, mb = mbar_w(iters)
res = summary(w, cv)
# PMFs
def pmf(x, w, bins):
    h, e = np.histogram(x, bins=bins, weights=w)
    with np.errstate(divide='ignore'):
        f = -kT * np.log(h / h.sum())
    f -= np.nanmin(f[np.isfinite(f)])
    return ((e[1:] + e[:-1]) / 2).tolist(), [float(v) if np.isfinite(v) else None for v in f]
res['pmf_Q_loop'] = pmf(cv[:, ix['Q_loop']], w, np.linspace(0, 1, 41))
res['pmf_phi43'] = pmf(cv[:, ix['phi43']], w, np.linspace(-180, 180, 37))
qx, qf = res['pmf_Q_loop']
qf_ = np.array([np.nan if v is None else v for v in qf]); qx_ = np.array(qx)
if np.isfinite(qf_[qx_ < 0.45]).any() and np.isfinite(qf_[qx_ >= 0.75]).any():
    mid = (qx_ >= 0.45) & (qx_ < 0.75)
    res['barrier_N_to_U_kcal'] = float(np.nanmax(qf_[mid]) - np.nanmin(qf_[qx_ >= 0.75])) if np.isfinite(qf_[mid]).any() else None
# block uncertainty
blocks = np.array_split(iters, 4)
bl = []
for b in blocks:
    if len(b) < 10:
        continue
    try:
        wb, cvb, _ = mbar_w(b); bl.append(summary(wb, cvb))
    except Exception as ex:
        print('block failed', ex)
def sd(key, sub=None):
    vals = [x[key][sub] if sub else x[key] for x in bl]
    vals = [v for v in vals if v is not None]
    return float(np.std(vals, ddof=1)) if len(vals) > 1 else None
res['block_sd'] = {'p_N': sd('p', 'N'), 'p_I': sd('p', 'I'), 'p_U': sd('p', 'U'), 'p_phi43_left': sd('p_phi43_left'),
                   'dG_I': sd('dG_kcal', 'I'), 'dG_U': sd('dG_kcal', 'U')}
# direct 310 K replica (no reweighting), as a cross-check
st0 = np.array([CV[t][np.where(S[t] == 0)[0][0]] for t in iters])
lab0, left0 = states(st0)
res['state0_direct'] = {'p': {s: float(np.mean(lab0 == s)) for s in ('N', 'I', 'U')}, 'p_phi43_left': float(np.mean(left0))}
acc = (z['n_acc'] / np.maximum(z['n_try'], 1)).tolist()
trips = 0
for r in range(K):
    last = None
    for t in range(n):
        s_ = S[t, r]
        if s_ == 0:
            if last == 'top': trips += 1
            last = 'bottom'
        elif s_ == K - 1 and last == 'bottom':
            last = 'top'
res.update(n_iter=int(n), burn_in=int(t0), g=g, n_samples=int(len(iters) * K), acceptance=acc, round_trips=trips,
           ns_per_replica=float(n * meta['steps_per_iter'] * meta['dt_ps'] / 1000), Teff=meta['Teff'], hot=meta['hot'])
json.dump(res, open(os.path.join(a.run, 'rest2_results.json'), 'w'), indent=1, default=float)
print(json.dumps({k: res[k] for k in ('p', 'dG_kcal', 'p_phi43_left', 'dG_phi43_left_vs_right_kcal', 'block_sd', 'round_trips', 'ns_per_replica')}, default=float))
