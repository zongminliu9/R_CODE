#!/usr/bin/env python3
"""Unbiased test: are the TUBA8-divergent positions enriched at any particular structural
class, or are they scattered? Permutation null = random position sets of the same size drawn
from the resolved alpha-tubulin residues, 100000 draws."""
import json, csv, numpy as np
ctx = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
rng = np.random.default_rng(7)

resolved = sorted(ctx)
divpos   = [d["canon_pos"] for d in div if d["canon_pos"] in ctx]
divuniq  = [d["canon_pos"] for d in div if d["canon_pos"] in ctx and d["tuba8_unique"]]
print(f"resolved alpha-tubulin residues in 9M1M: {len(resolved)} (auth {resolved[0]}-{resolved[-1]})")
print(f"TUBA8-divergent positions inside that range: {len(divpos)}  (TUBA8-unique: {len(divuniq)})")
print(f"divergent positions NOT resolved (C-term tail etc.): "
      f"{sorted(d['canon_pos'] for d in div if d['canon_pos'] not in ctx)}\n")

# structural classes
def cls(r):
    out = []
    if r["dSASA_cofactor"] > 5.0:  out.append("cofactor_iface")
    if r["dSASA_beta"]     > 5.0:  out.append("alphabeta_iface")
    if r["relSASA_mono"]   < 0.10: out.append("buried_core")
    if r["relSASA_mono"]   > 0.35: out.append("exposed_surface")
    if r.get("d_nucleotide", 99) < 5.0: out.append("nucleotide_pocket")
    if r["d_TBCE"] < 8.0: out.append("near_TBCE_8A")
    if r["d_TBCD"] < 8.0: out.append("near_TBCD_8A")
    return out

CLASSES = ["cofactor_iface","alphabeta_iface","buried_core","exposed_surface",
           "nucleotide_pocket","near_TBCE_8A","near_TBCD_8A"]
memb = {c: set() for c in CLASSES}
for p, r in ctx.items():
    for c in cls(r): memb[c].add(p)

NPERM = 100000
print(f"{'class':20s} {'|class|':>8s} {'obs':>5s} {'exp':>7s} {'ratio':>6s} {'p_enrich':>9s} {'p_deplete':>10s}")
print("-"*76)
res_arr = np.array(resolved)
out_rows = []
for c in CLASSES:
    S = memb[c]
    obs = sum(1 for p in divpos if p in S)
    # permutation
    hits = np.empty(NPERM, dtype=np.int32)
    mask = np.isin(res_arr, list(S))
    k = len(divpos)
    for i in range(NPERM):
        idx = rng.choice(len(res_arr), size=k, replace=False)
        hits[i] = mask[idx].sum()
    exp = hits.mean()
    p_up = (np.sum(hits >= obs)+1)/(NPERM+1)
    p_dn = (np.sum(hits <= obs)+1)/(NPERM+1)
    ratio = obs/exp if exp > 0 else float('nan')
    print(f"{c:20s} {len(S):8d} {obs:5d} {exp:7.2f} {ratio:6.2f} {p_up:9.5f} {p_dn:10.5f}")
    out_rows.append(dict(cls=c, n_class=len(S), observed=obs, expected=round(float(exp),3),
                         ratio=round(float(ratio),3), p_enrich=float(p_up), p_deplete=float(p_dn)))
json.dump(out_rows, open("results/structure/enrichment_permutation.json","w"), indent=1)

# ---- where IS the divergence concentrated? sliding-window density scan
print("\n=== divergence density scan (window = 15 residues, over FULL canonical 1-451) ===")
alldiv = np.array([d["canon_pos"] for d in div])
L = 451; W = 15
dens = np.array([((alldiv >= s) & (alldiv < s+W)).sum() for s in range(1, L-W+2)])
top = np.argsort(dens)[::-1][:8]
# window-level permutation: expected max density for 45 random positions
maxnull = np.empty(20000, dtype=np.int32)
for i in range(20000):
    pos = rng.choice(np.arange(1, L+1), size=len(alldiv), replace=False)
    maxnull[i] = max(((pos >= s) & (pos < s+W)).sum() for s in range(1, L-W+2, 3))
for t in sorted(top):
    s = t+1
    inwin = sorted(alldiv[(alldiv >= s) & (alldiv < s+W)].tolist())
    print(f"  residues {s:3d}-{s+W-1:3d}: {dens[t]:2d} substitutions  {inwin}")
print(f"\n  observed max window density = {dens.max()}")
print(f"  null distribution of max window density (45 random positions): "
      f"mean={maxnull.mean():.2f}, 99th pct={np.percentile(maxnull,99):.0f}, max={maxnull.max()}")
print(f"  p(max density >= {dens.max()}) = {(np.sum(maxnull >= dens.max())+1)/(len(maxnull)+1):.2e}")
