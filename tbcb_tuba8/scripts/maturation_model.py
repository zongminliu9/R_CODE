#!/usr/bin/env python3
"""Multistate alpha-tubulin maturation model with SHARED, LIMITING cofactor pools.

The point of sharing TBCB and TBCE between isoforms is that it makes the model
competitive: an isoform's fate depends not only on its own rate constants but on how
well it competes for a finite chaperone pool. That is what lets a single knockout
produce an isoform-selective AND region-selective phenotype, which no single-isoform
model can do.

Per isoform i:
    U_i  post-CCT, partly folded          (chaperone-committed)
    F_i  quasi-native, cofactor-free      (vulnerable)
    B_i  TBCB-bound
    E_i  TBCE-bound (post-handoff)
    H_i  mature soluble alpha/beta heterodimer
    M_i  assembled / retained polymer pool
    D_i  cleared (absorbing)

Shared:  TBCB_tot = TBCB_free + sum_i B_i ;  TBCE_tot = TBCE_free + sum_i E_i

Reactions (v) and stoichiometry (S) are built explicitly so the same code serves
simulation, steady-state solving, and sensitivity analysis.
"""
import numpy as np
from dataclasses import dataclass, field, asdict

STATES = ["U","F","B","E","H","M","D"]
NS = len(STATES)

@dataclass
class IsoformParams:
    name: str
    s:        float = 1.0     # synthesis flux into U (arb. units / s), scaled by region demand
    k_fold:   float = 1e-2    # U -> F   completion of folding after CCT release
    k_degU:   float = 1e-3    # U -> D   clearance of partly folded species
    kon_B:    float = 1.0     # F + TBCB -> B      (per uM per s)
    koff_B:   float = 1e-2    # B -> F + TBCB
    k_degF:   float = 1e-3    # F -> D   clearance of unprotected quasi-native species
    kon_Edir: float = 1e-2    # F + TBCE -> E      TBCB-INDEPENDENT capture (bypass route)
    k_tr:     float = 1.0     # B + TBCE -> E + TBCB   TBCB-mediated handoff
    k_degB:   float = 1e-5    # B -> D
    k_degE:   float = 1e-4    # E -> D
    k_asm:    float = 1.0     # E -> H + TBCE   heterodimer completion
    k_pol:    float = 1e-1    # H -> M
    k_depol:  float = 1e-2    # M -> H
    k_degH:   float = 1e-4    # H -> D
    k_degM:   float = 1e-5    # M -> D

@dataclass
class SystemParams:
    isoforms: list
    TBCB_tot: float = 1.0
    TBCE_tot: float = 1.0
    demand:   float = 1.0     # global multiplier on all synthesis fluxes (region tubulin demand)

def unpack(x, n_iso):
    """x = [states(n_iso*NS)] -> array (n_iso, NS)"""
    return x.reshape(n_iso, NS)

def free_cofactors(X, sp):
    B_bound = X[:, STATES.index("B")].sum()
    E_bound = X[:, STATES.index("E")].sum()
    TBCB_free = max(sp.TBCB_tot - B_bound, 0.0)
    TBCE_free = max(sp.TBCE_tot - E_bound, 0.0)
    return TBCB_free, TBCE_free

def rhs(t, x, sp):
    n = len(sp.isoforms)
    X = unpack(np.asarray(x, dtype=float), n)
    dX = np.zeros_like(X)
    iU, iF, iB, iE, iH, iM, iD = range(NS)
    Bf, Ef = free_cofactors(X, sp)
    for k, p in enumerate(sp.isoforms):
        U, F, B, E, H, M, D = X[k]
        U = max(U,0.0); F = max(F,0.0); B = max(B,0.0); E = max(E,0.0); H = max(H,0.0); M = max(M,0.0)
        v_syn   = p.s * sp.demand
        v_fold  = p.k_fold * U
        v_degU  = p.k_degU * U
        v_onB   = p.kon_B  * F * Bf
        v_offB  = p.koff_B * B
        v_degF  = p.k_degF * F
        v_onE   = p.kon_Edir * F * Ef
        v_tr    = p.k_tr   * B * Ef
        v_degB  = p.k_degB * B
        v_degE  = p.k_degE * E
        v_asm   = p.k_asm  * E
        v_pol   = p.k_pol  * H
        v_dep   = p.k_depol* M
        v_degH  = p.k_degH * H
        v_degM  = p.k_degM * M
        dX[k, iU] = v_syn - v_fold - v_degU
        dX[k, iF] = v_fold - v_onB + v_offB - v_degF - v_onE
        dX[k, iB] = v_onB - v_offB - v_tr - v_degB
        dX[k, iE] = v_onE + v_tr - v_asm - v_degE
        dX[k, iH] = v_asm - v_pol + v_dep - v_degH
        dX[k, iM] = v_pol - v_dep - v_degM
        dX[k, iD] = v_degU + v_degF + v_degB + v_degE + v_degH + v_degM
    return dX.reshape(-1)

def observable(X):
    """What an immunoblot of brain lysate reports: everything except cleared material and
    the pre-folded CCT-proximal species (which is a small, short-lived pool)."""
    iF, iB, iE, iH, iM = (STATES.index(s) for s in ("F","B","E","H","M"))
    return X[:, [iF, iB, iE, iH, iM]].sum(axis=1)

def steady_state(sp, t_end=5e6, rtol=1e-8, atol=1e-12):
    from scipy.integrate import solve_ivp
    n = len(sp.isoforms)
    x0 = np.zeros(n*NS)
    sol = solve_ivp(rhs, (0, t_end), x0, args=(sp,), method="LSODA",
                    rtol=rtol, atol=atol, dense_output=False)
    X = unpack(sol.y[:, -1], n)
    # verify stationarity of the non-absorbing states
    d = unpack(rhs(0.0, sol.y[:, -1], sp), n)
    drift = np.abs(d[:, :6]).max() / max(observable(X).max(), 1e-30)
    return X, sol, drift

if __name__ == "__main__":
    # smoke test: two isoforms, identical except TUBA8 lacks the TBCB-independent bypass
    a8  = IsoformParams("TUBA8",  s=0.15, kon_Edir=1e-3, k_degF=5e-3)
    a1a = IsoformParams("TUBA1A", s=1.00, kon_Edir=5e-2, k_degF=1e-3)
    for tbcb in (1.0, 0.0):
        sp = SystemParams([a8, a1a], TBCB_tot=tbcb, TBCE_tot=1.0, demand=1.0)
        X, sol, drift = steady_state(sp)
        obs = observable(X)
        print(f"TBCB_tot={tbcb}:  TUBA8={obs[0]:.4g}  TUBA1A={obs[1]:.4g}  (drift={drift:.1e})")
