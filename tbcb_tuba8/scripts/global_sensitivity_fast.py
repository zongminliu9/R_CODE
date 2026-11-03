#!/usr/bin/env python3
"""Q3 global sensitivity, using the exact fixed-point steady state.

Broad prior: each baseline rate constant over +/-1 decade, each TUBA8/TUBA1A ratio over
+/-1.5 decades, the three context variables (demand, TBCB_tot, TBCE_tot) over 3 decades, and
the TUBA8 synthesis flux over 2 decades below TUBA1A. 50000 draws.

Two outputs:
 (1) Spearman rho of selectivity against every sampled log-parameter -- which knobs matter
     across the whole prior, not at one hand-picked point.
 (2) The parameter region that REPRODUCES THE EXPERIMENT (TUBA8 strongly down, TUBA1A
     preserved), characterised by how each parameter's distribution shifts relative to the
     prior. This says what the observation implies, instead of asserting a mechanism.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from maturation_model import IsoformParams, SystemParams, observable, STATES
from exact_steady import steady_state_exact

rng = np.random.default_rng(20261103)
BASE = dict(k_fold=1e-2, k_degU=1e-3, kon_B=1.0, koff_B=1e-2, k_degF=1e-3,
            kon_Edir=1e-2, k_tr=1.0, k_degB=1e-5, k_degE=1e-4, k_asm=1.0,
            k_pol=1e-1, k_depol=1e-2, k_degH=1e-4, k_degM=1e-5)
RATIO_KEYS = ["k_fold","k_degU","kon_B","koff_B","k_degF","kon_Edir","k_tr",
              "k_degB","k_degE","k_asm","k_pol","k_degH"]
N = 50000

base_s = {k: BASE[k]*10**rng.uniform(-1, 1, N) for k in BASE}
lr     = {k: rng.uniform(-1.5, 1.5, N) for k in RATIO_KEYS}
demand = 10**rng.uniform(-1.5, 1.5, N)
TBCB   = 10**rng.uniform(-1.5, 1.5, N)
TBCE   = 10**rng.uniform(-1.5, 1.5, N)
ls8    = rng.uniform(-2.0, 0.0, N)

sel = np.full(N, np.nan); fc8 = np.full(N, np.nan); fc1 = np.full(N, np.nan)
bad = 0
for i in range(N):
    b  = {k: base_s[k][i] for k in BASE}
    p8 = dict(b); p1 = dict(b)
    for k in RATIO_KEYS: p8[k] = b[k]*10**lr[k][i]
    try:
        o = {}
        for tag, tb in (("WT", TBCB[i]), ("KO", 0.0)):
            sp = SystemParams([IsoformParams("TUBA8", s=10**ls8[i], **p8),
                               IsoformParams("TUBA1A", s=1.0, **p1)],
                              TBCB_tot=tb, TBCE_tot=TBCE[i], demand=demand[i])
            X, how, res = steady_state_exact(sp)
            if X is None or res > 1e-6: raise ValueError(how)
            o[tag] = observable(X)
        if min(o["WT"]) <= 0: raise ValueError("zero WT")
        f8 = o["KO"][0]/o["WT"][0]; f1 = o["KO"][1]/o["WT"][1]
        if not (np.isfinite(f8) and np.isfinite(f1) and f8 > 0 and f1 > 0): raise ValueError("nonfinite")
        fc8[i], fc1[i] = f8, f1
        sel[i] = np.log2(f8) - np.log2(f1)
    except Exception:
        bad += 1
        continue
    if i % 10000 == 0: print(f"  draw {i}/{N}", flush=True)

ok = np.isfinite(sel)
print(f"\nusable draws: {ok.sum()} / {N}  (rejected {bad})")
params = {}
for k in RATIO_KEYS: params["lr_"+k] = lr[k]
for k in BASE:       params["lb_"+k] = np.log10(base_s[k])
params["log_demand"] = np.log10(demand); params["log_TBCB"] = np.log10(TBCB)
params["log_TBCE"]   = np.log10(TBCE);   params["log_s8"]   = ls8

print("\n=== (1) Spearman rho of selectivity vs each sampled log-parameter ===")
print(f"{'parameter':22s} {'rho':>7s} {'p':>11s}")
rows = []
for k, v in params.items():
    rho, p = spearmanr(v[ok], sel[ok])
    rows.append((k, float(rho), float(p)))
rows.sort(key=lambda t: -abs(t[1]))
for k, rho, p in rows:
    star = "***" if p < 1e-20 else ("**" if p < 1e-6 else ("*" if p < 0.05 else ""))
    print(f"{k:22s} {rho:7.3f} {p:11.2e} {star}")

PHEN = ok & (sel < -1.5) & (fc8 < 0.5) & (fc1 > 0.75)
print(f"\n=== (2) draws reproducing the phenotype (sel<-1.5, FC8<0.5, FC1a>0.75): "
      f"{PHEN.sum()} / {ok.sum()} ({100*PHEN.sum()/ok.sum():.2f}%) ===")
print(f"{'parameter':22s} {'med|phen':>10s} {'med|rest':>10s} {'shift':>7s} {'p(MWU)':>11s}")
rows2 = []
rest = ok & ~PHEN
for k, v in params.items():
    a, b = v[PHEN], v[rest]
    _, p = mannwhitneyu(a, b)
    rows2.append((k, float(np.median(a)), float(np.median(b)),
                  float(np.median(a)-np.median(b)), float(p)))
rows2.sort(key=lambda t: t[4])
for k, ma, mb, sh, p in rows2:
    star = "***" if p < 1e-20 else ("**" if p < 1e-6 else ("*" if p < 0.05 else ""))
    print(f"{k:22s} {ma:10.3f} {mb:10.3f} {sh:+7.3f} {p:11.2e} {star}")

# how often is the phenotype reachable WITHOUT a kon_Edir deficit?
no_def = PHEN & (lr["kon_Edir"] > -0.3)
print(f"\n  phenotype draws with kon_Edir ratio > 0.5x (i.e. no real bypass deficit): "
      f"{no_def.sum()} / {PHEN.sum()} ({100*no_def.sum()/max(PHEN.sum(),1):.1f}%)")
for k in ["k_degF","kon_B","k_degE","koff_B"]:
    if no_def.sum() > 20:
        print(f"    among those, median lr_{k:9s} = {np.median(lr[k][no_def]):+.3f} "
              f"(prior median 0.000)")

os.makedirs("results/kinetics", exist_ok=True)
json.dump(dict(n_usable=int(ok.sum()), n_phenotype=int(PHEN.sum()), N=N,
               spearman=[dict(param=k, rho=r, p=p) for k, r, p in rows],
               phenotype_region=[dict(param=k, med_phen=a, med_rest=b, shift=s, p=p)
                                 for k, a, b, s, p in rows2],
               phenotype_without_bypass_deficit=int(no_def.sum())),
          open("results/kinetics/Q3_global_sensitivity.json","w"), indent=1)
np.savez_compressed("results/kinetics/Q3_raw.npz", sel=sel, fc8=fc8, fc1=fc1, **params)
print("\nwrote results/kinetics/Q3_global_sensitivity.json and Q3_raw.npz")
