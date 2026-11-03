#!/usr/bin/env python3
"""Fast steady state for the maturation network.

The cleared pool D is absorbing and grows without bound, so only the six live states have a
stationary point. We solve those directly with a damped Newton method on the residual, using
the conservation laws for the shared cofactor pools, and fall back to stiff integration if
Newton fails. This is ~100x faster than integrating to t=5e6 and makes large parameter sweeps
affordable."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maturation_model import STATES, NS, SystemParams, rhs, unpack
from scipy.optimize import fsolve

LIVE = [STATES.index(s) for s in ("U","F","B","E","H","M")]

def residual(y, sp, n):
    X = np.zeros((n, NS))
    X[:, LIVE] = np.maximum(y.reshape(n, len(LIVE)), 0.0)
    d = unpack(rhs(0.0, X.reshape(-1), sp), n)
    return d[:, LIVE].reshape(-1)

def steady_state_fast(sp, y0=None, verbose=False):
    n = len(sp.isoforms)
    if y0 is None:
        # crude but good initial guess: ignore competition, use linear-chain estimates
        g = []
        for p in sp.isoforms:
            U = p.s*sp.demand/(p.k_fold+p.k_degU)
            F = p.k_fold*U/max(p.kon_B*sp.TBCB_tot+p.k_degF+p.kon_Edir*sp.TBCE_tot, 1e-12)
            B = min(p.kon_B*F*sp.TBCB_tot/max(p.koff_B+p.k_tr*sp.TBCE_tot+p.k_degB,1e-12), sp.TBCB_tot)
            E = min((p.k_tr*B*sp.TBCE_tot+p.kon_Edir*F*sp.TBCE_tot)/max(p.k_asm+p.k_degE,1e-12), sp.TBCE_tot)
            H = p.k_asm*E/max(p.k_pol+p.k_degH,1e-12)
            M = p.k_pol*H/max(p.k_depol+p.k_degM,1e-12)
            g += [U,F,B,E,H,M]
        y0 = np.maximum(np.array(g), 1e-12)
    sol, info, ier, msg = fsolve(residual, y0, args=(sp, n), full_output=True,
                                 xtol=1e-12, maxfev=20000)
    y = np.maximum(sol, 0.0)
    r = np.abs(residual(y, sp, n))
    scale = max(y.max(), 1e-30)
    ok = (ier == 1) and (r.max()/scale < 1e-8)
    if not ok:
        from scipy.integrate import solve_ivp
        x0 = np.zeros(n*NS)
        s = solve_ivp(rhs, (0, 5e6), x0, args=(sp,), method="LSODA", rtol=1e-9, atol=1e-13)
        X = unpack(s.y[:, -1], n)
        return X, "integrated", float(np.abs(unpack(rhs(0,s.y[:,-1],sp),n)[:,LIVE]).max()/max(X[:,LIVE].max(),1e-30))
    X = np.zeros((n, NS)); X[:, LIVE] = y.reshape(n, len(LIVE))
    return X, "newton", float(r.max()/scale)

if __name__ == "__main__":
    import time
    from maturation_model import IsoformParams, steady_state, observable
    a8  = IsoformParams("TUBA8", s=0.15, kon_Edir=1e-3)
    a1a = IsoformParams("TUBA1A", s=1.0)
    sp = SystemParams([a8,a1a], TBCB_tot=1.0, TBCE_tot=1.0)
    t0=time.time(); Xi,_,di = steady_state(sp); t_int=time.time()-t0
    t0=time.time(); Xf,how,df = steady_state_fast(sp); t_new=time.time()-t0
    print(f"integration: {t_int*1000:8.1f} ms  obs={observable(Xi)}  drift={di:.1e}")
    print(f"newton     : {t_new*1000:8.1f} ms  obs={observable(Xf)}  ({how}, res={df:.1e})")
    rel = np.abs(observable(Xf)-observable(Xi))/np.maximum(np.abs(observable(Xi)),1e-30)
    print(f"max relative disagreement on observable: {rel.max():.3e}   speedup {t_int/t_new:.0f}x")
