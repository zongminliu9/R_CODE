#!/usr/bin/env python3
"""Which physically plausible rate differences make TBCB loss hit TUBA8 selectively,
and why would that be confined to a high-demand region such as cerebellum?

Four questions, each answered by a computation rather than an argument:

 Q0 NULL MODEL     Can low TUBA8 expression alone, with NO isoform-specific rate constants,
                   produce selective loss purely through mass-action competition?
 Q1 DECOMPOSITION  Taking one candidate mechanism at a time, how large must the TUBA8/TUBA1A
                   ratio be before the knockout becomes selective?
 Q2 REGIONALITY    Scanning synthesis demand against cofactor abundance, is there a regime
                   where loss is severe at high demand and absent at low demand?
 Q3 SENSITIVITY    Over a broad Latin-hypercube prior, which parameters actually control
                   selectivity, and which are irrelevant?

Selectivity is reported as  sel = log2(FC_TUBA8) - log2(FC_TUBA1A),  FC = KO/WT observable.
A strongly negative sel means TUBA8 drops while TUBA1A does not.
"""
import json, os, csv, itertools, warnings, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from maturation_model import IsoformParams, SystemParams, steady_state, observable, STATES

os.makedirs("results/kinetics", exist_ok=True)
rng = np.random.default_rng(20251103)

BASE = dict(k_fold=1e-2, k_degU=1e-3, kon_B=1.0, koff_B=1e-2, k_degF=1e-3,
            kon_Edir=1e-2, k_tr=1.0, k_degB=1e-5, k_degE=1e-4, k_asm=1.0,
            k_pol=1e-1, k_depol=1e-2, k_degH=1e-4, k_degM=1e-5)

def run(a8_over=None, a1_over=None, s8=0.15, s1=1.0, TBCB=1.0, TBCE=1.0, demand=1.0):
    p8 = dict(BASE); p8.update(a8_over or {})
    p1 = dict(BASE); p1.update(a1_over or {})
    out = {}
    for tag, tb in (("WT", TBCB), ("KO", 0.0)):
        sp = SystemParams([IsoformParams("TUBA8", s=s8, **p8),
                           IsoformParams("TUBA1A", s=s1, **p1)],
                          TBCB_tot=tb, TBCE_tot=TBCE, demand=demand)
        X, _, drift = steady_state(sp)
        out[tag] = (observable(X), X, drift)
    o_wt, o_ko = out["WT"][0], out["KO"][0]
    eps = 1e-300
    fc8 = (o_ko[0]+eps)/(o_wt[0]+eps); fc1 = (o_ko[1]+eps)/(o_wt[1]+eps)
    sel = np.log2(fc8) - np.log2(fc1)
    return dict(fc_TUBA8=fc8, fc_TUBA1A=fc1, sel=sel,
                wt8=o_wt[0], wt1=o_wt[1], ko8=o_ko[0], ko1=o_ko[1],
                drift=max(out["WT"][2], out["KO"][2]))

# ---------------------------------------------------------------- Q0 null model
print("="*78)
print("Q0  NULL MODEL -- identical rate constants, TUBA8 differs ONLY in synthesis flux")
print("="*78)
print(f"{'s_TUBA8/s_TUBA1A':>18s} {'FC_TUBA8':>10s} {'FC_TUBA1A':>10s} {'sel (log2)':>11s}")
q0 = []
for ratio in [1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
    r = run(s8=ratio, s1=1.0)
    q0.append(dict(ratio=ratio, **r))
    print(f"{ratio:18.3g} {r['fc_TUBA8']:10.4f} {r['fc_TUBA1A']:10.4f} {r['sel']:11.4f}")
json.dump(q0, open("results/kinetics/Q0_null_model.json","w"), indent=1, default=float)

# ---------------------------------------------------------------- Q1 decomposition
print()
print("="*78)
print("Q1  ONE MECHANISM AT A TIME -- how big a TUBA8/TUBA1A ratio is needed?")
print("="*78)
MECHS = {
 "kon_Edir (TBCB-independent TBCE capture)": ("kon_Edir", [1.0, 0.5, 0.2, 0.1, 0.03, 0.01, 0.003]),
 "k_degF   (clearance of unprotected F)":    ("k_degF",   [1.0, 2.0, 5.0, 10., 30., 100., 300.]),
 "k_fold   (completion of folding)":         ("k_fold",   [1.0, 0.5, 0.2, 0.1, 0.03, 0.01, 0.003]),
 "kon_B    (TBCB capture rate)":             ("kon_B",    [1.0, 2.0, 5.0, 10., 30., 100., 300.]),
 "k_tr     (TBCB->TBCE handoff)":            ("k_tr",     [1.0, 0.5, 0.2, 0.1, 0.03, 0.01, 0.003]),
 "k_asm    (heterodimer completion)":        ("k_asm",    [1.0, 0.5, 0.2, 0.1, 0.03, 0.01, 0.003]),
}
q1 = {}
for label, (key, ratios) in MECHS.items():
    print(f"\n  {label}")
    print(f"    {'ratio':>8s} {'FC_TUBA8':>10s} {'FC_TUBA1A':>10s} {'sel':>9s}")
    rec = []
    for rt in ratios:
        r = run(a8_over={key: BASE[key]*rt})
        rec.append(dict(ratio=rt, **r))
        print(f"    {rt:8.3g} {r['fc_TUBA8']:10.4f} {r['fc_TUBA1A']:10.4f} {r['sel']:9.3f}")
    q1[label] = rec
json.dump(q1, open("results/kinetics/Q1_mechanism_decomposition.json","w"), indent=1, default=float)
