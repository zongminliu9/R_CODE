"""Umbrella-sampling PMF along one backbone phi dihedral (default residue 43) in the folded alpha monomer, with
Hamiltonian exchange between neighbouring windows (periodic: window 0 and window K-1 are also neighbours through
the ladder ends only).  Bias 0.5*k*d^2 with d the periodic angle difference, k = 200 kJ/mol/rad^2, 24 windows every
15 degrees.  Analysis: MBAR -> W(phi) in kcal/mol, dG(phi>0 basin vs phi<0 basin) and the barriers between them.
This is the free-energy cost for residue 43 to adopt the left-handed (Gly-like) backbone seen at 43 in the
cofactor-bound alpha structures, in each sequence context.
Usage: a17_us_phi.py run <build> <run_dir> <rep> [--res 43] [--ns 4]   |   a17_us_phi.py analyse <run_dir>/r<k>
"""
import os, sys, json, argparse, math
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('mode'); ap.add_argument('a1'); ap.add_argument('a2', nargs='?'); ap.add_argument('rep', nargs='?', type=int, default=0)
ap.add_argument('--res', type=int, default=43); ap.add_argument('--ns', type=float, default=4.0)
ap.add_argument('--k', type=float, default=200.0); ap.add_argument('--nwin', type=int, default=24)
ap.add_argument('--budget', type=float, default=2800.0); ap.add_argument('--platform', default=None)
ap.add_argument('--quick', action='store_true', help='smoke test: tiny equilibration and starts')
ap.add_argument('--steps', type=int, default=250)
a = ap.parse_args()
kT = 0.0019872043 * 310.0
if a.mode == 'run':
    import openmm as mm
    from openmm import app, unit
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import common as C
    from hrex import HREX
    W = os.path.join(a.a2, f'r{a.rep}'); os.makedirs(W, exist_ok=True)
    LOCAL = os.path.join(os.environ.get('L_SCRATCH_JOB') or os.environ.get('TMPDIR') or '/tmp', f'usphi_{os.path.basename(a.a2)}_r{a.rep}')
    os.makedirs(LOCAL, exist_ok=True)
    log = lambda m: print(m, flush=True)
    top, pos = C.load_cif(os.path.join(a.a1, 'minimized.cif.gz'))
    byres = {(int(x.residue.id), x.name): x.index for x in top.atoms() if x.residue.name not in C.SOLVENT | {'GTP', 'MG'}}
    q = [byres[(a.res - 1, 'C')], byres[(a.res, 'N')], byres[(a.res, 'CA')], byres[(a.res, 'C')]]
    seed = 6060 + 19 * a.rep
    EQ = os.path.join(W, 'equil_state.npz')
    if not os.path.exists(EQ):
        P, V, box, rec = C.equilibrate(top, pos, T=310.0, seed=seed, platform=a.platform, log=log, npt=True, scale=0.01 if a.quick else 1.0)
        json.dump(rec, open(os.path.join(W, 'equil.json'), 'w'), indent=1)
        C.save_state_npz(EQ, P, V, box)
    P, V, box, _ = C.load_state_npz(EQ)
    system = C.make_system(top, hmr=True)
    b = mm.CustomTorsionForce('0.5*kphi*dd^2; dd = dt - 2*pi*floor((dt+pi)/(2*pi)); dt = theta - phi0; pi = 3.141592653589793')
    b.addGlobalParameter('kphi', a.k); b.addGlobalParameter('phi0', 0.0)
    b.addTorsion(*q, [])
    system.addForce(b)
    centres = np.linspace(-math.pi, math.pi, a.nwin, endpoint=False)
    states = [{'phi0': float(c), 'kphi': a.k} for c in centres]
    def phi(X):
        return math.radians(C.dihedral(*[X[i] for i in q]))
    # starting structures: drive phi through the windows from the equilibrated value (slow rotation, 0.2 ns/window)
    STARTS = os.path.join(W, 'starts.npz')
    if not os.path.exists(STARTS):
        integ = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        plat, props = C.get_platform(a.platform)
        ctx = mm.Context(system, integ, plat, props)
        ctx.setPeriodicBoxVectors(*box); ctx.setPositions(P); ctx.setVelocities(V)
        p0 = phi(P); order = np.argsort(np.abs(((centres - p0) + math.pi) % (2 * math.pi) - math.pi))
        starts = {}
        cur = p0
        # walk outward in both directions from the window nearest the current phi
        for direction in (+1, -1):
            ctx.setPositions(P); ctx.setVelocities(V); cur = p0
            k0 = int(order[0])
            for step in range(a.nwin // 2 + 1):
                k = (k0 + direction * step) % a.nwin
                tgt = centres[k]
                dlt = ((tgt - cur) + math.pi) % (2 * math.pi) - math.pi
                for x in np.linspace(cur, cur + dlt, 21)[1:]:
                    ctx.setParameter('phi0', float(((x + math.pi) % (2 * math.pi)) - math.pi)); integ.step(5 if a.quick else 250)
                cur = cur + dlt
                integ.step(50 if a.quick else 25000)
                X = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(C.NM)
                if k not in starts:
                    starts[k] = X.copy(); log(f'window {k}: target {math.degrees(tgt):.0f}, phi {math.degrees(phi(X)):.1f}')
        np.savez(STARTS, pos=np.array([starts[k] for k in range(a.nwin)]))
        del ctx, integ            # release the GPU before the HREX context (see a15)
        import gc; gc.collect()
    z = np.load(STARTS)
    n_iter = int(round(a.ns * 1000 / (a.steps * 0.004)))
    h = HREX(system, z['pos'], box, states, temperature_K=310.0, timestep_ps=0.004, steps_per_iter=a.steps,
             platform=a.platform, shared_dir=W, local_dir=LOCAL, tag='usphi', seed=seed, cv_fn=lambda X: [phi(X)])
    h.resume_or_init()
    done = h.run(n_iter, a.budget, ck_every_min=10.0, log_every=100)
    json.dump({'iteration': h.iteration, 'target': n_iter, 'done': bool(done), 'acceptance': h.acceptance(),
               'centres': centres.tolist(), 'res': a.res, 'k': a.k}, open(os.path.join(W, 'usphi_status.json'), 'w'), indent=1)
else:
    import pymbar
    W = a.a1
    z = np.load(os.path.join(W, 'usphi_ck.npz'), allow_pickle=True); stt = json.load(open(os.path.join(W, 'usphi_status.json')))
    U = np.array(z['rec_u']); S = np.array(z['rec_s']); CV = np.array([np.array(x, float) for x in z['rec_cv']])[:, :, 0]
    n, K, _ = U.shape; t0 = int(0.1 * n); its = np.arange(t0, n)
    u_kn = np.concatenate([U[t].T for t in its], axis=1); st = np.concatenate([S[t] for t in its]); ph = np.concatenate([CV[t] for t in its])
    N_k = np.bincount(st, minlength=K); centres = np.array(stt['centres']); kk = stt['k'] / (0.0083144626 * 310.0)
    d = (ph - centres[st] + math.pi) % (2 * math.pi) - math.pi
    u0 = u_kn[st, np.arange(len(st))] - 0.5 * kk * d ** 2
    fes = pymbar.FES(u_kn, N_k, mbar_options={'solver_protocol': 'robust'})
    edges = np.linspace(-math.pi, math.pi, 37)
    fes.generate_fes(u0, ph, fes_type='histogram', histogram_parameters={'bin_edges': edges})
    mid = 0.5 * (edges[1:] + edges[:-1])
    r = fes.get_fes(mid, reference_point='from-lowest', uncertainty_method='analytical')
    Wp = np.array(r['f_i']) * kT; dW = np.array(r['df_i']) * kT
    w = np.exp(-Wp / kT); w = np.where(np.isfinite(w), w, 0)
    pos_ = mid > 0
    dG = -kT * math.log(w[pos_].sum() / w[~pos_].sum())
    out = dict(phi_deg=np.degrees(mid).tolist(), W_kcal=[float(x) if np.isfinite(x) else None for x in Wp],
               dW_kcal=[float(x) if np.isfinite(x) else None for x in dW], dG_left_minus_right_kcal=float(dG),
               p_left=float(w[pos_].sum() / w.sum()), res=stt['res'], n_iter=int(n),
               acceptance=(z['n_acc'] / np.maximum(z['n_try'], 1)).tolist())
    # barriers on the two paths between the basins (through phi = 0 and through phi = 180)
    def seg_max(lo, hi):
        m = (mid >= math.radians(lo)) & (mid <= math.radians(hi)); return float(np.nanmax(Wp[m]))
    out['barrier_via_0_kcal'] = seg_max(-40, 40) - float(np.nanmin(Wp))
    out['barrier_via_180_kcal'] = max(seg_max(140, 180), seg_max(-180, -140)) - float(np.nanmin(Wp))
    json.dump(out, open(os.path.join(W, 'usphi_results.json'), 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('dG_left_minus_right_kcal', 'p_left', 'barrier_via_0_kcal', 'barrier_via_180_kcal')}))
