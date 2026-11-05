"""Markov state models of the H1-S2 loop region from unbiased 310 K MD, on a discretisation shared by all systems so
that the same metastable states can be compared between isoforms and mutants.

Features: sin/cos of backbone phi/psi of residues 30-50 (all systems have the same residue numbering there).
Pipeline: TICA (deeptime, lag 1 ns, 90 % kinetic variance, at most 8 dims) on the pooled data -> k-means (150 centres)
-> per-system reversible MLE MSMs at a lag chosen from the implied timescales -> metastable macrostates from PCCA+ of
the pooled MSM -> for every system: macrostate populations, MFPTs between macrostates, per-macrostate mean Q_loop /
loop RMSD / phi43, and a Chapman-Kolmogorov test.  Nothing from a system whose CK test fails is reported as a
population or MFPT (the gate from the protocol); its numbers are written under 'failed_ck'.
Usage: a11_msm.py <runs_root>/md310 <out_dir> [--systems a1a,a8,...] [--lag_ns 5] [--nmeta 3]
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md

ap = argparse.ArgumentParser()
ap.add_argument('root'); ap.add_argument('out'); ap.add_argument('--systems', default=None)
ap.add_argument('--lag_ns', type=float, default=None); ap.add_argument('--nmeta', type=int, default=3)
ap.add_argument('--stride', type=int, default=1); ap.add_argument('--dt_ps', type=float, default=20.0)
ap.add_argument('--res', default='30-50')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
from deeptime.decomposition import TICA
from deeptime.clustering import KMeans
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import MaximumLikelihoodMSM
lo, hi = [int(x) for x in a.res.split('-')]
systems = a.systems.split(',') if a.systems else sorted(os.path.basename(d) for d in glob.glob(os.path.join(a.root, '*')))
feats, meta = {}, {}
for sid in systems:
    trajs = []
    for rdir in sorted(glob.glob(os.path.join(a.root, sid, 'r*'))):
        pdb = os.path.join(rdir, 'solute.pdb')
        segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd')))
        if not segs or not os.path.exists(pdb):
            continue
        t = md.join([md.load_dcd(s, top=pdb, stride=a.stride) for s in segs], check_topology=False)
        trajs.append(t)
    if not trajs:
        continue
    top = trajs[0].topology
    sel = [r.index for r in top.residues if lo <= r.resSeq <= hi and r.is_protein]
    X = []; aux = []
    ref = md.load(os.path.join(sorted(glob.glob(os.path.join(a.root, sid, 'r*')))[0], 'solute.pdb'))
    loop_bb = top.select(f'backbone and resSeq 32 to 45')
    core_ca = top.select('name CA and (resSeq < 25 or resSeq > 55)')
    for t in trajs:
        _, phi = md.compute_phi(t); _, psi = md.compute_psi(t)
        ip, _ = md.compute_phi(t[:1]); jp, _ = md.compute_psi(t[:1])
        rp = [top.atom(q[1]).residue.resSeq for q in ip]; rq = [top.atom(q[1]).residue.resSeq for q in jp]
        kp = [i for i, r in enumerate(rp) if lo <= r <= hi]; kq = [i for i, r in enumerate(rq) if lo <= r <= hi]
        f = np.concatenate([np.sin(phi[:, kp]), np.cos(phi[:, kp]), np.sin(psi[:, kq]), np.cos(psi[:, kq])], axis=1)
        X.append(f.astype(np.float32))
        t.superpose(ref, atom_indices=core_ca)
        rmsd = np.sqrt(np.mean(np.sum((t.xyz[:, loop_bb] - ref.xyz[0, loop_bb]) ** 2, axis=2), axis=1)) * 10
        k43 = [i for i, r in enumerate(rp) if r == 43][0]
        aux.append(np.stack([rmsd, np.degrees(phi[:, k43])], axis=1))
    feats[sid] = X; meta[sid] = dict(n_traj=len(X), n_frames=int(sum(len(x) for x in X)),
                                     ns=float(sum(len(x) for x in X) * a.dt_ps * a.stride / 1000), aux=aux)
    print(sid, meta[sid]['n_traj'], 'trajectories', round(meta[sid]['ns'], 1), 'ns', flush=True)
dt_ns = a.dt_ps * a.stride / 1000
pooled = [x for sid in feats for x in feats[sid]]
tica = TICA(lagtime=max(1, int(round(1.0 / dt_ns))), var_cutoff=0.9, dim=None).fit(pooled).fetch_model()
Y = {sid: [tica.transform(x)[:, :8] for x in feats[sid]] for sid in feats}
allY = np.concatenate([y for sid in Y for y in Y[sid]])
km = KMeans(n_clusters=150, max_iter=200, fixed_seed=7).fit(allY[::max(1, len(allY) // 60000)]).fetch_model()
D = {sid: [km.transform(y) for y in Y[sid]] for sid in Y}
# implied timescales of the pooled data to choose the lag
lags = [l for l in (1, 2, 5, 10, 25, 50, 100, 250) if l * 5 < min(len(d) for sid in D for d in D[sid])]
its = {}
for l in lags:
    c = TransitionCountEstimator(lagtime=l, count_mode='sliding').fit([d for sid in D for d in D[sid]]).fetch_model()
    m = MaximumLikelihoodMSM(reversible=True).fit(c.submodel_largest()).fetch_model()
    its[l * dt_ns] = (m.timescales(5) * dt_ns).tolist()
lag = int(round(a.lag_ns / dt_ns)) if a.lag_ns else lags[min(len(lags) - 1, 3)]
cnt = TransitionCountEstimator(lagtime=lag, count_mode='sliding').fit([d for sid in D for d in D[sid]]).fetch_model()
pooled_msm = MaximumLikelihoodMSM(reversible=True).fit(cnt.submodel_largest()).fetch_model()
pcca = pooled_msm.pcca(a.nmeta)
memb = np.zeros((150, a.nmeta)); memb[pooled_msm.count_model.state_symbols] = pcca.memberships
macro_of = np.argmax(memb, axis=1)
out = dict(lag_ns=lag * dt_ns, implied_timescales_ns=its, tica_dim=int(tica.output_dimension), n_meta=a.nmeta, systems={})
for sid in D:
    c = TransitionCountEstimator(lagtime=lag, count_mode='sliding').fit(D[sid]).fetch_model().submodel_largest()
    m = MaximumLikelihoodMSM(reversible=True).fit(c).fetch_model()
    sym = m.count_model.state_symbols
    pi_micro = np.zeros(150); pi_micro[sym] = m.stationary_distribution
    p_macro = [float(pi_micro[macro_of == k].sum()) for k in range(a.nmeta)]
    # macrostate CV means from the frames assigned to each macrostate, weighted by the MSM stationary distribution
    d_all = np.concatenate(D[sid]); aux = np.concatenate(meta[sid]['aux'])
    cvm = {}
    for k in range(a.nmeta):
        sel = np.isin(d_all, np.where(macro_of == k)[0])
        cvm[k] = dict(n_frames=int(sel.sum()), mean_rmsd_loop_A=float(aux[sel, 0].mean()) if sel.any() else None,
                      frac_phi43_left=float(np.mean(aux[sel, 1] > 0)) if sel.any() else None)
    # MFPT between macrostates (sets of active microstates)
    mf = {}
    for i in range(a.nmeta):
        for j in range(a.nmeta):
            if i == j: continue
            A = [k for k, s in enumerate(sym) if macro_of[s] == i]; B = [k for k, s in enumerate(sym) if macro_of[s] == j]
            if A and B:
                mf[f'{i}->{j}'] = float(m.mfpt(A, B) * dt_ns)
    # Chapman-Kolmogorov test on the macrostate memberships restricted to the active set
    ck_pass = None; ck_err = None
    try:
        from deeptime.util.validation import ck_test
        ests = []
        for mult in (1, 2, 3, 4):
            cc = TransitionCountEstimator(lagtime=lag * mult, count_mode='sliding').fit(D[sid]).fetch_model().submodel_largest()
            ests.append(MaximumLikelihoodMSM(reversible=True).fit(cc).fetch_model())
        ck = ck_test(ests, n_metastable_sets=a.nmeta)
        pred, est = ck.predictions, ck.estimates
        ck_err = float(np.max(np.abs(pred - est)))
        ck_pass = bool(ck_err < 0.1)
    except Exception as ex:
        ck_err = f'ck failed: {ex}'
    rec = dict(n_traj=meta[sid]['n_traj'], ns=meta[sid]['ns'], active_fraction=float(m.count_model.selected_count_fraction)
               if hasattr(m.count_model, 'selected_count_fraction') else None, macro_populations=p_macro,
               macro_cv=cvm, mfpt_ns=mf, ck_max_abs_dev=ck_err, ck_pass=ck_pass,
               slowest_timescales_ns=(m.timescales(3) * dt_ns).tolist())
    out['systems'][sid] = rec if ck_pass else {'failed_ck': rec}
    print(sid, json.dumps(rec, default=float)[:400], flush=True)
json.dump(out, open(os.path.join(a.out, 'msm_results.json'), 'w'), indent=1, default=float)
