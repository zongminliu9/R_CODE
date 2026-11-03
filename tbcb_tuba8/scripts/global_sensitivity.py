#!/usr/bin/env python3
"""Q3 global sensitivity.

Rather than trusting one hand-picked parameter set, we draw a broad prior: every baseline rate
constant varies over +/-1 decade, every TUBA8/TUBA1A ratio over +/-1.5 decades, and the three
context variables (demand, TBCB_tot, TBCE_tot) over 3 decades. Then we ask two things:

  (1) Which log-ratios control selectivity across the whole prior (Spearman rho)?
  (2) Among draws that REPRODUCE THE EXPERIMENT -- TUBA8 strongly down, TUBA1A preserved --
      how is each parameter distributed relative to the prior? That characterises the
      parameter region the observation actually implies, instead of asserting a mechanism.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
from maturation_model import IsoformParams, SystemParams, steady_state, observable, STATES

rng = np.random.default_rng(424242)
BASE = dict(k_fold=1e-2, k_degU=1e-3, kon_B=1.0, koff_B=1e-2, k_degF=1e-3,
            kon_Edir=1e-2, k_tr=1.0, k_degB=1e-5, k_degE=1e-4, k_asm=1.0,
            k_pol=1e-1, k_depol=1e-2, k_degH=1e-4, k_degM=1e-5)
RATIO_KEYS = ["k_fold","k_degU","kon_B","koff_B","k_degF","kon_Edir","k_tr",
              "k_degB","k_degE","k_asm","k_pol","k_degH"]

N = 1500
recs = []
for it in range(N):
    base = {k: v*10**rng.uniform(-1, 1) for k, v in BASE.items()}
    lr   = {k: rng.uniform(-1.5, 1.5) for k in RATIO_KEYS}
    demand = 10**rng.uniform(-1.5, 1.5)
    TBCB   = 10**rng.uniform(-1.5, 1.5)
    TBCE   = 10**rng.uniform(-1.5, 1.5)
    ls8    = rng.uniform(-2.0, 0.0)          # log10(s_TUBA8 / s_TUBA1A)
    p8 = dict(base); p1 = dict(base)
    for k in RATIO_KEYS: p8[k] = base[k]*10**lr[k]
    try:
        out = {}
        for tag, tb in (("WT", TBCB), ("KO", 0.0)):
            sp = SystemParams([IsoformParams("TUBA8", s=10**ls8, **p8),
                               IsoformParams("TUBA1A", s=1.0, **p1)],
                              TBCB_tot=tb, TBCE_tot=TBCE, demand=demand)
            X, _, drift = steady_state(sp, t_end=2e7)
            out[tag] = (observable(X), drift)
        e = 1e-300
        o_wt, o_ko = out["WT"][0], out["KO"][0]
        if min(o_wt) < 1e-20: continue
        fc8 = (o_ko[0]+e)/(o_wt[0]+e); fc1 = (o_ko[1]+e)/(o_wt[1]+e)
        if not (np.isfinite(fc8) and np.isfinite(fc1)): continue
        rec = dict(sel=float(np.log2(fc8)-np.log2(fc1)), fc8=float(fc8), fc1=float(fc1),
                   log_demand=np.log10(demand), log_TBCB=np.log10(TBCB), log_TBCE=np.log10(TBCE),
                   log_s8=ls8, drift=float(max(out["WT"][1], out["KO"][1])))
        for k in RATIO_KEYS: rec["lr_"+k] = lr[k]
        for k in BASE: rec["lb_"+k] = float(np.log10(base[k]))
        recs.append(rec)
    except Exception:
        continue

print(f"usable draws: {len(recs)} / {N}")
keys = [k for k in recs[0] if k.startswith(("lr_","lb_","log_"))]
sel  = np.array([r["sel"] for r in recs])
fc8  = np.array([r["fc8"] for r in recs]); fc1 = np.array([r["fc1"] for r in recs])

print("\n=== (1) Spearman rho of selectivity vs each sampled log-parameter ===")
rows=[]
for k in keys:
    v = np.array([r[k] for r in recs])
    rho, p = spearmanr(v, sel)
    rows.append((k, float(rho), float(p)))
rows.sort(key=lambda t: -abs(t[1]))
print(f"{'parameter':22s} {'rho':>7s} {'p':>11s}")
for k, rho, p in rows:
    star = "***" if p < 1e-10 else ("**" if p < 1e-4 else ("*" if p < 0.05 else ""))
    print(f"{k:22s} {rho:7.3f} {p:11.2e} {star}")

# ---- (2) parameter region consistent with the experiment
PHEN = (sel < -1.5) & (fc8 < 0.5) & (fc1 > 0.75)
print(f"\n=== (2) draws reproducing the phenotype (sel<-1.5, FC_TUBA8<0.5, FC_TUBA1A>0.75): "
      f"{PHEN.sum()} / {len(recs)}  ({100*PHEN.mean():.1f}%) ===")
print(f"{'parameter':22s} {'median|phen':>12s} {'median|rest':>12s} {'shift':>7s} {'p(MWU)':>10s}")
rows2=[]
for k in keys:
    v = np.array([r[k] for r in recs])
    if PHEN.sum() < 8: break
    a, b = v[PHEN], v[~PHEN]
    try: _, p = mannwhitneyu(a, b)
    except Exception: p = 1.0
    rows2.append((k, float(np.median(a)), float(np.median(b)), float(np.median(a)-np.median(b)), float(p)))
rows2.sort(key=lambda t: t[4])
for k, ma, mb, sh, p in rows2:
    star = "***" if p < 1e-10 else ("**" if p < 1e-4 else ("*" if p < 0.05 else ""))
    print(f"{k:22s} {ma:12.3f} {mb:12.3f} {sh:7.3f} {p:10.2e} {star}")

json.dump(dict(n=len(recs), n_phen=int(PHEN.sum()),
               spearman=[dict(param=k, rho=r, p=p) for k,r,p in rows],
               phenotype_region=[dict(param=k, med_phen=a, med_rest=b, shift=s, p=p)
                                 for k,a,b,s,p in rows2]),
          open("results/kinetics/Q3_global_sensitivity.json","w"), indent=1)
