#!/usr/bin/env python3
"""Q2 regionality and Q3 global sensitivity for the maturation model.

Q2 asks the question the experiment forces on us: the knockout phenotype is confined to
cerebellum, yet TUBA8's amino-acid sequence is the same in every region. So the regional
term must be a FLUX or ABUNDANCE variable, not a residue. We scan synthesis demand against
cofactor abundance and look for a regime that is severe at high demand and silent at low.

Q3 draws parameters from a broad log-uniform prior and asks which ones actually control
selectivity, so the conclusion does not rest on one hand-picked parameter set.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from maturation_model import IsoformParams, SystemParams, steady_state, observable, STATES

os.makedirs("results/kinetics", exist_ok=True)
rng = np.random.default_rng(31337)
BASE = dict(k_fold=1e-2, k_degU=1e-3, kon_B=1.0, koff_B=1e-2, k_degF=1e-3,
            kon_Edir=1e-2, k_tr=1.0, k_degB=1e-5, k_degE=1e-4, k_asm=1.0,
            k_pol=1e-1, k_depol=1e-2, k_degH=1e-4, k_degM=1e-5)

def evaluate(p8, p1, s8, s1, TBCB, TBCE, demand):
    res = {}
    for tag, tb in (("WT", TBCB), ("KO", 0.0)):
        sp = SystemParams([IsoformParams("TUBA8", s=s8, **p8), IsoformParams("TUBA1A", s=s1, **p1)],
                          TBCB_tot=tb, TBCE_tot=TBCE, demand=demand)
        X, _, drift = steady_state(sp)
        iF,iB,iE,iH,iM = (STATES.index(x) for x in ("F","B","E","H","M"))
        res[tag] = dict(total=observable(X), sol=X[:, [iF,iB,iE,iH]].sum(axis=1),
                        pol=X[:, iM], drift=drift)
    e = 1e-300
    fc  = (res["KO"]["total"]+e)/(res["WT"]["total"]+e)
    fcs = (res["KO"]["sol"]+e)/(res["WT"]["sol"]+e)
    return dict(fc8=float(fc[0]), fc1=float(fc[1]), sel=float(np.log2(fc[0])-np.log2(fc[1])),
                fc8_sol=float(fcs[0]), fc1_sol=float(fcs[1]),
                sel_sol=float(np.log2(fcs[0])-np.log2(fcs[1])),
                wt8=float(res["WT"]["total"][0]), wt1=float(res["WT"]["total"][1]),
                drift=float(max(res["WT"]["drift"], res["KO"]["drift"])))

# ---------------------------------------------------------------- Q2 regionality
print("="*96)
print("Q2  REGIONALITY: does a high-demand region become vulnerable while a low-demand one stays silent?")
print("    TUBA8 carries ONE difference: kon_Edir reduced 30x (deficient TBCB-independent route).")
print("="*96)
p8 = dict(BASE); p8["kon_Edir"] = BASE["kon_Edir"]/30.0
p1 = dict(BASE)
demands = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
print(f"{'demand':>8s} {'FC_TUBA8':>9s} {'FC_TUBA1A':>10s} {'sel':>7s} | "
      f"{'FC8_sol':>8s} {'FC1_sol':>8s} {'sel_sol':>8s} | {'WT8':>10s}")
q2 = []
for d in demands:
    r = evaluate(p8, p1, 0.15, 1.0, 1.0, 1.0, d)
    q2.append(dict(demand=d, **r))
    print(f"{d:8.3g} {r['fc8']:9.4f} {r['fc1']:10.4f} {r['sel']:7.3f} | "
          f"{r['fc8_sol']:8.4f} {r['fc1_sol']:8.4f} {r['sel_sol']:8.3f} | {r['wt8']:10.4g}")
json.dump(q2, open("results/kinetics/Q2_regionality_demand.json","w"), indent=1)

print("\n  --- 2D scan: demand x TBCE abundance (sel on total protein) ---")
tbces = [0.1, 0.3, 1.0, 3.0, 10.0]
print(f"{'':>10s}" + "".join(f"{'TBCE='+str(t):>11s}" for t in tbces))
grid = []
for d in demands:
    row = []
    for t in tbces:
        r = evaluate(p8, p1, 0.15, 1.0, 1.0, t, d)
        row.append(r["sel"]); grid.append(dict(demand=d, TBCE=t, **r))
    print(f"demand={d:<4.3g}" + "".join(f"{v:11.3f}" for v in row))
json.dump(grid, open("results/kinetics/Q2_grid_demand_TBCE.json","w"), indent=1)
