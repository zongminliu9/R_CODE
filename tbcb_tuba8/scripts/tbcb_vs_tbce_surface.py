#!/usr/bin/env python3
"""Does AF3 put TBCB on the same alpha-tubulin surface that TBCE occupies experimentally?

This is a consistency check on the AF3 model that uses an independent experimental structure
rather than AF3's own confidence metrics, which are calibrated on PDB-like end states and give
no guarantee for a chaperone bound to a non-native substrate.

It also bears on mechanism. If TBCB and TBCE compete for the same alpha-tubulin surface, the
TBCB->TBCE handoff must be a displacement, and the kinetic model's transfer step k_tr acquires a
physical meaning. If the surfaces are disjoint, both could bind at once and the handoff is
something else entirely.
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np

af3 = json.load(open("results/af3/interface_aligned.json"))
ctx = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
div = {d["canon_pos"] for d in json.load(open("results/alignment/tuba8_divergent_positions.json"))}

tbcb = set(af3["TUBA8"]["robust_canon"]) | set(af3["TUBA1A"]["robust_canon"])
# experimental partners in 9M1M, by buried surface and proximity
tbce = {p for p,r in ctx.items() if r["dSASA_cofactor"] > 5 and r["d_TBCE"] < 5.0}
tbcd = {p for p,r in ctx.items() if r["dSASA_cofactor"] > 5 and r["d_TBCD"] < 5.0}
beta = {p for p,r in ctx.items() if r["dSASA_beta"] > 5}
cof  = {p for p,r in ctx.items() if r["dSASA_cofactor"] > 5}
resolved = set(ctx)

print(f"alpha-tubulin residues resolved in 9M1M: {len(resolved)} (1-437)")
print(f"AF3 TBCB interface: {len(tbcb)} residues; "
      f"{len(tbcb & resolved)} of them fall inside the resolved range")
print(f"experimental 9M1M interfaces: TBCE {len(tbce)}, TBCD {len(tbcd)}, beta-tubulin {len(beta)}")

def jac(a, b):
    a, b = set(a) & resolved, set(b) & resolved
    return len(a & b)/max(len(a | b), 1), len(a & b)

print(f"\n{'pair':34s} {'overlap':>8s} {'Jaccard':>8s}")
for nm, s in (("TBCB(AF3) vs TBCE(9M1M)", tbce), ("TBCB(AF3) vs TBCD(9M1M)", tbcd),
              ("TBCB(AF3) vs beta-tubulin", beta), ("TBCB(AF3) vs any cofactor", cof)):
    j, n = jac(tbcb, s)
    print(f"{nm:34s} {n:8d} {j:8.3f}")

shared = sorted((tbcb & tbce))
print(f"\nresidues shared by the AF3 TBCB interface and the experimental TBCE interface: "
      f"{shared if shared else 'NONE'}")
print(f"TBCB-only (within resolved range): {sorted((tbcb & resolved) - cof)[:40]}")

# ---- is the AF3 TBCB interface on a surface that is FREE in 9M1M?
free = (tbcb & resolved) - cof - beta
print(f"\nAF3 TBCB interface residues that are unoccupied in 9M1M: {len(free)}/"
      f"{len(tbcb & resolved)}  ({100*len(free)/max(len(tbcb&resolved),1):.0f}%)")

# ---- domain decomposition: where is TBCB binding, and where is TUBA8 divergent?
DOMAINS = {"N-terminal GTPase (1-205)": range(1,206),
           "intermediate (206-381)":    range(206,382),
           "C-terminal (382-440)":      range(382,441),
           "C-terminal tail (441-451)":  range(441,452)}
print(f"\n{'domain':28s} {'size':>5s} {'TBCB iface':>11s} {'frac':>6s} "
      f"{'TUBA8 div':>10s} {'frac':>6s}")
rows=[]
for nm, rng in DOMAINS.items():
    R = set(rng); n = len(R)
    ti = len(tbcb & R); dv = len(div & R)
    print(f"{nm:28s} {n:5d} {ti:11d} {ti/n:6.3f} {dv:10d} {dv/n:6.3f}")
    rows.append(dict(domain=nm, size=n, tbcb_interface=ti, tbcb_frac=ti/n,
                     divergent=dv, div_frac=dv/n))
json.dump(dict(tbcb_iface=sorted(tbcb), tbce_iface=sorted(tbce), tbcd_iface=sorted(tbcd),
               overlap_tbcb_tbce=shared, free_in_9M1M=sorted(free), domains=rows),
          open("results/af3/tbcb_vs_tbce_surface.json","w"), indent=1)
