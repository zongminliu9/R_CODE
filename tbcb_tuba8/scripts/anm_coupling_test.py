#!/usr/bin/env python3
"""Two questions the ANM can answer rigorously from a single structure.

 (1) Do the TUBA8-divergent positions occupy mechanically special sites -- unusually high
     effectiveness (they push others), sensitivity (they get pushed) or mobility?
 (2) The 32-45 loop physically touches no cofactor. Does it nevertheless couple
     mechanically to the functional surfaces -- the TBCE contact patch, the nucleotide
     pocket, the alpha-beta interface, the C-terminal region?

For (2) the null is every other contiguous 14-residue window in the chain, which controls for
window length and for the general tendency of exposed loops to be mobile.
"""
import json, numpy as np
from scipy.stats import mannwhitneyu

prs  = np.load("results/anm/prs_matrix.npy")
anm  = json.load(open("results/anm/alpha_anm_prs.json"))
ctx  = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
div  = json.load(open("results/alignment/tuba8_divergent_positions.json"))
pos  = np.array([r["canon_pos"] for r in anm]); N = len(pos)
p2i  = {p: i for i, p in enumerate(pos)}
rng  = np.random.default_rng(99)

divpos = [d["canon_pos"] for d in div if d["canon_pos"] in p2i]
loop   = [p for p in range(32, 46) if p in p2i]

print(f"ANM residues {N}; divergent inside {len(divpos)}; 32-45 loop residues present {len(loop)}")

# ---------------- (1) are divergent positions mechanically special?
print("\n=== (1) mechanical character of the divergent set (permutation, 1e5 draws) ===")
print(f"{'metric':16s} {'mean|div':>10s} {'mean|rest':>10s} {'p_high':>8s} {'p_low':>8s}")
NP = 100000
idx_div = np.array([p2i[p] for p in divpos])
for key in ("eff_z", "sen_z", "msf_z"):
    v = np.array([r[key] for r in anm])
    obs = v[idx_div].mean()
    rest = np.delete(v, idx_div).mean()
    null = np.array([v[rng.choice(N, size=len(idx_div), replace=False)].mean() for _ in range(NP)])
    p_hi = (np.sum(null >= obs)+1)/(NP+1); p_lo = (np.sum(null <= obs)+1)/(NP+1)
    print(f"{key:16s} {obs:10.4f} {rest:10.4f} {p_hi:8.5f} {p_lo:8.5f}")

# ---------------- (2) directed coupling from the 32-45 loop to functional surfaces
TARGETS = {
 "TBCE contact (dSASA>5)":  [p for p,r in ctx.items() if r["dSASA_cofactor"]>5 and r["d_TBCE"]<5.0],
 "any cofactor contact":    [p for p,r in ctx.items() if r["dSASA_cofactor"]>5],
 "nucleotide pocket <5A":   [p for p,r in ctx.items() if r.get("d_nucleotide",99)<5.0],
 "alpha-beta interface":    [p for p,r in ctx.items() if r["dSASA_beta"]>5],
 "C-terminal 420-437":      [p for p in range(420,438) if p in p2i],
}
print("\n=== (2) mechanical coupling of the 32-45 loop to surfaces it does not touch ===")
print("    null = all other contiguous 14-residue windows (same length)")
print(f"{'target surface':26s} {'|T|':>4s} {'coupling':>10s} {'null mean':>10s} {'pct':>6s} {'p':>8s}")
W = len(loop)
windows = [list(range(s, s+W)) for s in range(int(pos.min()), int(pos.max())-W+2)]
windows = [w for w in windows if all(p in p2i for p in w)]
rows=[]
for tname, T in TARGETS.items():
    T = [p for p in T if p in p2i]
    if not T: continue
    ti = np.array([p2i[p] for p in T])
    def coup(wlist):
        wi = np.array([p2i[p] for p in wlist])
        sub = prs[np.ix_(wi, ti)]
        # exclude self/neighbour pairs: coupling to residues within 4 of the window is trivial
        mask = np.abs(pos[wi][:,None] - pos[ti][None,:]) > 4
        return sub[mask].mean() if mask.any() else np.nan
    obs = coup(loop)
    nullv = np.array([coup(w) for w in windows])
    nullv = nullv[np.isfinite(nullv)]
    pct = 100.0*np.mean(nullv < obs)
    p   = (np.sum(nullv >= obs)+1)/(len(nullv)+1)
    print(f"{tname:26s} {len(T):4d} {obs:10.4g} {nullv.mean():10.4g} {pct:6.1f} {p:8.4f}")
    rows.append(dict(target=tname, n_target=len(T), coupling=float(obs),
                     null_mean=float(nullv.mean()), percentile=float(pct), p=float(p)))
json.dump(rows, open("results/anm/loop32_45_coupling.json","w"), indent=1)

# ---------------- which single residues are the strongest effectors onto the TBCE patch?
T = [p for p,r in ctx.items() if r["dSASA_cofactor"]>5 and p in p2i]
ti = np.array([p2i[p] for p in T])
infl = np.array([prs[p2i[p], ti][np.abs(pos[ti]-p)>4].mean() for p in pos])
order = np.argsort(infl)[::-1]
print(f"\n=== strongest single-residue effectors onto the cofactor-contact patch (|T|={len(T)}) ===")
print("    rank  pos  influence   divergent?  region")
reg = {d["canon_pos"]: d["region"] for d in div}
for k in range(15):
    p = int(pos[order[k]])
    print(f"    {k+1:4d} {p:4d} {infl[order[k]]:10.4g}   "
          f"{'YES' if p in reg else '-':>4s}        {reg.get(p,'')}")
dz = (infl-infl.mean())/infl.std()
print(f"\n  mean influence z-score of divergent positions: "
      f"{dz[[p2i[p] for p in divpos]].mean():+.4f}")
