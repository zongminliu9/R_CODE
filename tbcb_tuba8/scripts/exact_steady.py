#!/usr/bin/env python3
"""Exact, fast steady state by fixed-point iteration on the free-cofactor concentrations.

The network is bilinear only through the three terms kon_B*F*Bf, kon_Edir*F*Ef and k_tr*B*Ef.
Hold the free pools (Bf, Ef) fixed and each isoform's chain becomes LINEAR, with a closed-form
stationary solution. Then update (Bf, Ef) from the conservation laws and iterate. This replaces
stiff integration to t=5e6 and is orders of magnitude faster, which is what makes a broad
parameter sweep affordable.

Closed form at fixed (Bf, Ef), per isoform:
  U = s*demand / (k_fold + k_degU)
  B = kon_B*Bf*F / (koff_B + k_tr*Ef + k_degB)                        == betaB * F
  F = k_fold*U / ( kon_B*Bf + k_degF + kon_Edir*Ef - koff_B*betaB )
  E = Ef*(kon_Edir*F + k_tr*B) / (k_asm + k_degE)
  H = k_asm*E / ( k_pol + k_degH - k_depol*k_pol/(k_depol + k_degM) )
  M = k_pol*H / (k_depol + k_degM)
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from maturation_model import STATES, NS

def _chain(p, demand, Bf, Ef):
    U = p.s*demand/(p.k_fold + p.k_degU)
    denB = p.koff_B + p.k_tr*Ef + p.k_degB
    betaB = p.kon_B*Bf/denB if denB > 0 else 0.0
    denF = p.kon_B*Bf + p.k_degF + p.kon_Edir*Ef - p.koff_B*betaB
    if denF <= 1e-300:          # degenerate: no loss channel for F
        return None
    F = p.k_fold*U/denF
    B = betaB*F
    E = Ef*(p.kon_Edir*F + p.k_tr*B)/(p.k_asm + p.k_degE)
    denH = p.k_pol + p.k_degH - p.k_depol*p.k_pol/(p.k_depol + p.k_degM)
    if denH <= 1e-300:
        return None
    H = p.k_asm*E/denH
    M = p.k_pol*H/(p.k_depol + p.k_degM)
    return np.array([U, F, B, E, H, M, 0.0])

def steady_state_exact(sp, tol=1e-13, maxit=4000):
    n = len(sp.isoforms)
    Bf, Ef = sp.TBCB_tot, sp.TBCE_tot
    X = np.zeros((n, NS))
    for it in range(maxit):
        rows = [_chain(p, sp.demand, Bf, Ef) for p in sp.isoforms]
        if any(r is None for r in rows): return None, "degenerate", np.inf
        Xn = np.array(rows)
        Bsum = Xn[:, STATES.index("B")].sum()
        Esum = Xn[:, STATES.index("E")].sum()
        # the bound amount cannot exceed the pool; solve the scaling self-consistently
        Bf_new = sp.TBCB_tot/(1.0 + Bsum/max(Bf, 1e-300)) if Bf > 0 else 0.0
        Ef_new = sp.TBCE_tot/(1.0 + Esum/max(Ef, 1e-300)) if Ef > 0 else 0.0
        w = 0.5                                    # damping keeps the iteration stable
        Bf2 = (1-w)*Bf + w*Bf_new
        Ef2 = (1-w)*Ef + w*Ef_new
        err = max(abs(Bf2-Bf)/max(sp.TBCB_tot, 1e-300) if sp.TBCB_tot > 0 else 0.0,
                  abs(Ef2-Ef)/max(sp.TBCE_tot, 1e-300))
        Bf, Ef = Bf2, Ef2
        X = Xn
        if err < tol: break
    # verify against the conservation laws we claimed to satisfy
    Bsum = X[:, STATES.index("B")].sum(); Esum = X[:, STATES.index("E")].sum()
    resid = max(abs(sp.TBCB_tot - (Bf + Bsum))/max(sp.TBCB_tot, 1e-30) if sp.TBCB_tot > 0
                else abs(Bsum),
                abs(sp.TBCE_tot - (Ef + Esum))/max(sp.TBCE_tot, 1e-30))
    return X, f"fixedpoint({it+1})", float(resid)

if __name__ == "__main__":
    import time
    from maturation_model import IsoformParams, SystemParams, steady_state, observable
    print(f"{'case':34s} {'ODE obs':>24s} {'exact obs':>24s} {'rel.diff':>10s} {'speedup':>8s}")
    cases = [
      ("baseline", dict(kon_Edir=1e-3), 1.0, 1.0, 1.0),
      ("TBCB KO", dict(kon_Edir=1e-3), 0.0, 1.0, 1.0),
      ("high demand", dict(kon_Edir=1e-3), 1.0, 1.0, 30.0),
      ("scarce TBCE", dict(kon_Edir=1e-3), 1.0, 0.1, 1.0),
      ("abundant TBCE", dict(kon_Edir=1e-2), 1.0, 10.0, 3.0),
    ]
    for lab, ov, tb, te, dem in cases:
        a8 = IsoformParams("TUBA8", s=0.15, **ov); a1 = IsoformParams("TUBA1A", s=1.0)
        sp = SystemParams([a8, a1], TBCB_tot=tb, TBCE_tot=te, demand=dem)
        t0=time.time(); Xi,_,_ = steady_state(sp); t1=time.time()-t0
        t0=time.time(); Xe,how,res = steady_state_exact(sp); t2=time.time()-t0
        oi, oe = observable(Xi), observable(Xe)
        rd = np.max(np.abs(oe-oi)/np.maximum(np.abs(oi),1e-30))
        print(f"{lab:34s} {np.array2string(oi,precision=4,max_line_width=200):>24s} "
              f"{np.array2string(oe,precision=4,max_line_width=200):>24s} {rd:10.2e} {t1/t2:7.0f}x")
