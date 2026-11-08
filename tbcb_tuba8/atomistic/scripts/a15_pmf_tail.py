"""Absolute binding free energy of the alpha tail to a CAP-Gly groove from an umbrella-sampling PMF (with Hamiltonian
exchange between neighbouring windows), giving dG0 and Kd with an explicit standard-state term.

CV r = distance between the CA centroid of the groove residues and the N/CA/C centroid of the tail's C-terminal
residue.  Windows every 0.1 nm from r_bound - 0.1 to r_bound + 2.4 nm, k = 1000 kJ/mol/nm^2; window starting
structures from a slow pull (0.5 nm/ns).  MBAR (a15 analyse) gives W(r); with no orientational restraint the bulk
satisfies W(r) = -kT ln(4 pi r^2) + C, and
    K = integral over the site of exp(-(W(r) - C)/kT) dr,   dG0 = -kT ln(K C0),  C0 = 1/1660.54 A^3,
with the site boundary at the barrier maximum between bound and bulk.  The bulk constant C is fitted on the last
0.6 nm of windows; the fit residual is reported (if W does not follow -2kT ln r there, the box did not reach bulk and
dG0 is not reported).
Usage: a15_pmf_tail.py run <big_box_build> <run_dir> <rep> --groove A:187,188,203,204,223 --term B:451 [--ns 8]
       a15_pmf_tail.py analyse <run_dir>/r<k>
"""
import os, sys, json, time, argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('mode'); ap.add_argument('a1'); ap.add_argument('a2', nargs='?'); ap.add_argument('rep', nargs='?', type=int, default=0)
ap.add_argument('--groove', default='A:187,188,203,204,223'); ap.add_argument('--term', default='B:451')
ap.add_argument('--ns', type=float, default=8.0); ap.add_argument('--budget', type=float, default=2800.0)
ap.add_argument('--k', type=float, default=1000.0); ap.add_argument('--span', type=float, default=2.4)
ap.add_argument('--spacing', type=float, default=0.1); ap.add_argument('--platform', default=None)
ap.add_argument('--quick', action='store_true', help='smoke test: tiny equilibration and starts')
ap.add_argument('--steps', type=int, default=500)
ap.add_argument('--g1names', default='CA'); ap.add_argument('--g2names', default='N,CA,C')
ap.add_argument('--start_offset', type=float, default=-0.1, help='first window relative to the starting distance (nm)')
a = ap.parse_args()
kT_kj = 0.0083144626 * 310.0
if a.mode == 'run':
    import openmm as mm
    from openmm import app, unit
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import common as C
    from hrex import HREX
    W = os.path.join(a.a2, f'r{a.rep}'); os.makedirs(W, exist_ok=True)
    LOCAL = os.path.join(os.environ.get('L_SCRATCH_JOB') or os.environ.get('TMPDIR') or '/tmp', f'pmf_{os.path.basename(a.a2)}_r{a.rep}')
    os.makedirs(LOCAL, exist_ok=True)
    log = lambda m: print(m, flush=True)
    top, pos = C.load_cif(os.path.join(a.a1, 'minimized.cif.gz'))
    def sel(spec, names):
        cid, ids = spec.split(':'); out = set()
        for tok in ids.split(','):
            if '-' in tok:
                lo_, hi_ = tok.split('-'); out |= set(range(int(lo_), int(hi_) + 1))
            else:
                out.add(int(tok))
        return [x.index for x in top.atoms() if (cid == '*' or x.residue.chain.id == cid) and x.residue.name not in C.SOLVENT
                and int(x.residue.id) in out and x.name in names]
    g1 = sel(a.groove, tuple(a.g1names.split(','))); g2 = sel(a.term, tuple(a.g2names.split(',')))
    assert g1 and g2
    seed = 5150 + 13 * a.rep
    EQ = os.path.join(W, 'equil_state.npz')
    if not os.path.exists(EQ):
        P, V, box, rec = C.equilibrate(top, pos, T=310.0, seed=seed, platform=a.platform, log=log, npt=True, scale=0.01 if a.quick else 1.0)
        json.dump(rec, open(os.path.join(W, 'equil.json'), 'w'), indent=1)
        C.save_state_npz(EQ, P, V, box)
    P, V, box, _ = C.load_state_npz(EQ)
    def cv(X):
        return float(np.linalg.norm(X[g1].mean(0) - X[g2].mean(0)))
    system = C.make_system(top, hmr=True)
    [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0].setUseDispersionCorrection(False)
    us = mm.CustomCentroidBondForce(2, '0.5*kus*(distance(g1,g2)-r0)^2')
    us.addGlobalParameter('kus', a.k); us.addGlobalParameter('r0', cv(P))
    us.addGroup(g1); us.addGroup(g2); us.addBond([0, 1], []); us.setUsesPeriodicBoundaryConditions(True)
    system.addForce(us)
    rb = cv(P)
    centres = np.round(np.arange(rb + a.start_offset, rb + a.span + 1e-9, a.spacing), 4)
    PULL = os.path.join(W, 'pull_starts.npz')
    if not os.path.exists(PULL):
        integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds)
        plat, props = C.get_platform(a.platform)
        ctx = mm.Context(system, integ, plat, props)
        ctx.setPeriodicBoxVectors(*box); ctx.setPositions(P); ctx.setVelocities(V)
        starts = {}
        rr = rb + a.start_offset
        ctx.setParameter('r0', rr); integ.step(50 if a.quick else 5000)
        nstep = int(round(a.spacing / 0.5 * 1000 / 0.004))      # 0.5 nm/ns
        for c in centres:
            while rr < c - 1e-6:
                rr = min(c, rr + 0.002); ctx.setParameter('r0', rr); integ.step(1 if a.quick else max(1, nstep // 50))
            integ.step(25 if a.quick else 2500)
            X = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(C.NM)
            starts[float(c)] = (X.copy(), cv(X))
            log(f'pull: window {c:.2f} nm, r = {cv(X):.3f}')
        np.savez(PULL, centres=centres, pos=np.array([starts[float(c)][0] for c in centres]),
                 r=np.array([starts[float(c)][1] for c in centres]))
        # release the GPU: on Sherlock a second CUDA context in this process fails ("No compatible CUDA device")
        del ctx, integ
        import gc; gc.collect()
    z = np.load(PULL)
    states = [{'r0': float(c), 'kus': a.k} for c in z['centres']]
    n_iter = int(round(a.ns * 1000 / (a.steps * 0.004)))
    h = HREX(system, z['pos'], box, states, temperature_K=310.0, timestep_ps=0.004, steps_per_iter=a.steps,
             platform=a.platform, shared_dir=W, local_dir=LOCAL, tag='pmf', seed=seed, cv_fn=lambda X: [cv(X)])
    h.resume_or_init()
    done = h.run(n_iter, a.budget, ck_every_min=10.0, log_every=100)
    json.dump({'iteration': h.iteration, 'target': n_iter, 'done': bool(done), 'acceptance': h.acceptance(),
               'centres': z['centres'].tolist(), 'r_bound_start': rb, 'groove': a.groove, 'term': a.term},
              open(os.path.join(W, 'pmf_status.json'), 'w'), indent=1)
    log(f'iteration {h.iteration}/{n_iter} done={done}')
else:
    import pymbar
    W = a.a1
    z = np.load(os.path.join(W, 'pmf_ck.npz'), allow_pickle=True)
    stt = json.load(open(os.path.join(W, 'pmf_status.json')))
    U = np.array(z['rec_u']); S = np.array(z['rec_s']); CVr = np.array([np.array(x, float) for x in z['rec_cv']])[:, :, 0]
    n, K, _ = U.shape
    t0 = int(0.1 * n)
    its = np.arange(t0, n)
    u_kn = np.concatenate([U[t].T for t in its], axis=1); st = np.concatenate([S[t] for t in its]); r = np.concatenate([CVr[t] for t in its])
    N_k = np.bincount(st, minlength=K)
    centres = np.array(stt['centres'])
    # unbiased reduced potential = u_k - bias_k for the sampled state of each sample
    bias = 0.5 * a.k * (r - centres[st]) ** 2 / kT_kj
    u0 = u_kn[st, np.arange(len(st))] - bias
    mb = pymbar.MBAR(u_kn, N_k, solver_protocol='robust')
    fes = pymbar.FES(u_kn, N_k, mbar_options={'solver_protocol': 'robust'})
    edges = np.linspace(centres.min() - 0.1, centres.max() + 0.1, 2 * K + 3)
    fes.generate_fes(u0, r, fes_type='histogram', histogram_parameters={'bin_edges': edges})
    mid = 0.5 * (edges[1:] + edges[:-1])
    # pymbar's histogram FES has no entry for an empty bin (KeyError (0,) on 2026-11-07), so it is evaluated only in
    # bins that hold samples; empty bins stay NaN and are counted (an empty bin inside the range means poor overlap)
    has = np.histogram(r, bins=edges)[0] > 0
    res = fes.get_fes(mid[has], reference_point='from-lowest', uncertainty_method='analytical')
    Wr = np.full(len(mid), np.nan); dW = np.full(len(mid), np.nan)
    Wr[has] = np.array(res['f_i']) * kT_kj / 4.184; dW[has] = np.array(res['df_i']) * kT_kj / 4.184   # kcal/mol
    n_empty = int((~has).sum()); n_out = int(np.sum((r < edges[0]) | (r > edges[-1])))
    kT = kT_kj / 4.184
    ok = np.isfinite(Wr)
    bulk = ok & (mid > mid[ok].max() - 0.6)
    # bulk model W = -kT ln(4 pi r^2) + C
    Cfit = np.mean(Wr[bulk] + kT * np.log(4 * np.pi * mid[bulk] ** 2))
    resid = float(np.std(Wr[bulk] + kT * np.log(4 * np.pi * mid[bulk] ** 2)))
    imin = int(np.nanargmin(np.where(ok, Wr, np.nan)))
    # site boundary: maximum of W between the minimum and the bulk region
    after = np.where(ok & (mid > mid[imin]) & (mid < mid[bulk].min()))[0]
    ib = int(after[np.argmax(Wr[after])]) if len(after) else imin + 3
    site = ok & (mid <= mid[ib])
    dr_nm = mid[1] - mid[0]
    K_nm3 = float(np.sum(np.exp(-(Wr[site] - Cfit) / kT)) * dr_nm)           # nm^3
    K_A3 = K_nm3 * 1000.0
    dG0 = -kT * np.log(K_A3 / 1660.54)
    Kd_M = 1.0 / (K_A3 / 1660.54)
    out = dict(r_nm=mid.tolist(), W_kcal=[float(x) if np.isfinite(x) else None for x in Wr], dW_kcal=[float(x) if np.isfinite(x) else None for x in dW],
               bulk_fit_C_kcal=float(Cfit), bulk_fit_residual_kcal=resid, r_min_nm=float(mid[imin]), r_site_boundary_nm=float(mid[ib]),
               barrier_kcal=float(Wr[ib] - Wr[imin]), K_A3=K_A3, dG0_kcal=float(dG0), Kd_M=float(Kd_M),
               gate_bulk_reached=bool(resid < 0.3), acceptance=(z['n_acc'] / np.maximum(z['n_try'], 1)).tolist(), n_iter=int(n),
               n_empty_bins=n_empty, n_samples_outside_edges=n_out,
               note='dG0 and Kd are only valid if gate_bulk_reached and the window overlap is adequate')
    ov = mb.compute_overlap()['matrix']; out['min_adjacent_overlap'] = float(min(ov[k, k + 1] for k in range(K - 1)))
    # ---- derived kinetics: position-dependent diffusion from window residence segments (Hummer 2005,
    # D = var(r)/tau_ac in a harmonic window), then the 1-D Smoluchowski MFPT out of the site (reflecting at the
    # inner edge, absorbing 0.3 nm beyond the site boundary); koff = 1/MFPT, kon = koff/Kd.  1-D projection: an
    # order-of-magnitude estimate, labelled 'derived'.
    dt_ps = a.steps * 0.004
    Dw = {}
    for k in range(K):
        segs = []
        for rep in range(K):
            cur = []
            for t in range(t0, n):
                if S[t, rep] == k:
                    cur.append(CVr[t, rep])
                else:
                    if len(cur) >= 20: segs.append(np.array(cur))
                    cur = []
            if len(cur) >= 20: segs.append(np.array(cur))
        if not segs:
            continue
        num = 0.0; den = 0.0
        for x in segs:
            x = x - x.mean(); v = np.mean(x * x)
            if v <= 0: continue
            ac = [1.0]
            for lag in range(1, min(len(x) // 2, 50)):
                c = np.mean(x[:-lag] * x[lag:]) / v
                if c < 0.05: break
                ac.append(c)
            tau = (0.5 + sum(ac[1:])) * dt_ps
            num += v / tau * len(x); den += len(x)
        if den:
            Dw[float(centres[k])] = num / den          # nm^2/ps
    if len(Dw) >= 3:
        cs = np.array(sorted(Dw)); Ds = np.array([Dw[c] for c in cs])
        Dr = np.interp(mid, cs, Ds)
        r_out = mid[ib] + 0.3
        sel = ok & (mid <= r_out)
        Wsel = Wr[sel] / kT; rr = mid[sel]; Dsel = Dr[sel]
        inner = np.cumsum(np.exp(-Wsel)) * dr_nm
        integrand = np.exp(Wsel) / Dsel * inner
        start = np.searchsorted(rr, mid[imin])
        mfpt_ps = float(np.sum(integrand[start:]) * dr_nm)
        koff = 1.0 / (mfpt_ps * 1e-12)
        out.update(D_nm2_ps={str(k): v for k, v in Dw.items()}, mfpt_unbind_ps=mfpt_ps, koff_s=koff,
                   kon_M_s=koff / Kd_M, residence_time_s=mfpt_ps * 1e-12,
                   kinetics_note='derived: 1-D Smoluchowski on the computed PMF with window-estimated D(r)')
    json.dump(out, open(os.path.join(W, 'pmf_results.json'), 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('dG0_kcal', 'Kd_M', 'barrier_kcal', 'bulk_fit_residual_kcal', 'gate_bulk_reached', 'min_adjacent_overlap')}))
