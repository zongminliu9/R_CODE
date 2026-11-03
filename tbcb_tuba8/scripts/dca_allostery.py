#!/usr/bin/env python3
"""The allostery question, asked properly.

A strong direct coupling between two residues that TOUCH is not allostery -- it is a contact.
The rigorous signature of functional long-range coupling is a pair that is strongly coupled in
evolution while being FAR APART in space, because no local physical interaction explains it.

So we restrict to pairs with Cb-Cb > 12 A in 9M1M and ask:
  (a) is the 32-45 loop coupled, through such non-contacting pairs, to the functional surfaces
      (the AF3-predicted TBCB interface, the nucleotide pocket, the alpha-beta interface)?
  (b) do TUBA8's 45 divergent positions coevolve with each other more than chance, i.e. do they
      behave as a module rather than as independent drift?
  (c) does the coevolutionary answer agree with the elastic-network answer, or not?

Null model throughout: all other contiguous windows of the same length, which controls for
window size and for each region's intrinsic coupling density.
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx
from scipy.stats import mannwhitneyu

DI = np.load("results/dca/dca_APC_full.npy").astype(np.float64)
L = DI.shape[0]
meta = json.load(open("results/dca/dca_meta_full.json"))
print(f"DCA: L={L}, M={meta['M']}, Meff={meta['Meff']:.1f} (Meff/L={meta['Meff_over_L']:.2f})")

# structure distances (canonical numbering), Cb-Cb
arr = pdbx.get_structure(pdbx.CIFFile.read("data/pdb/9M1M.cif"), model=1, use_author_fields=True)
A = arr[(arr.chain_id=="a") & struc.filter_amino_acids(arr)]
pos = {}
for rid in np.unique(A.res_id):
    s = A[A.res_id==rid]
    cb = s[s.atom_name=="CB"]; ca = s[s.atom_name=="CA"]
    p = cb.coord[0] if cb.array_length() else (ca.coord[0] if ca.array_length() else None)
    if p is not None: pos[int(rid)] = p
D = np.full((L, L), np.nan)
ks = sorted(pos)
for a_ in ks:
    for b_ in ks:
        D[a_-1, b_-1] = np.linalg.norm(pos[a_]-pos[b_])
resolved = np.array([k-1 for k in ks])

ctx = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
af3 = json.load(open("results/af3/interface_aligned.json"))
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
DIV = sorted({d["canon_pos"] for d in div})

TARGETS = {
 "TBCB interface (AF3)":   sorted(set(af3["TUBA8"]["robust_canon"]) | set(af3["TUBA1A"]["robust_canon"])),
 "late-cofactor interface":[p for p,r in ctx.items() if r["dSASA_cofactor"] > 5],
 "nucleotide pocket":      [p for p,r in ctx.items() if r.get("d_nucleotide",99) < 5],
 "alpha-beta interface":   [p for p,r in ctx.items() if r["dSASA_beta"] > 5],
}
LOOP = list(range(32, 46))
FAR = 12.0

def noncontact_coupling(src, tgt):
    """mean APC-corrected coupling over pairs that are >FAR apart in space (and >=6 in sequence)."""
    vals = []
    for i in src:
        for j in tgt:
            if abs(i-j) < 6: continue
            d = D[i-1, j-1]
            if not np.isfinite(d) or d <= FAR: continue
            vals.append(DI[i-1, j-1])
    return (np.mean(vals) if vals else np.nan), len(vals)

print(f"\n=== (a) is the 32-45 loop coupled to functional surfaces it does not touch? ===")
print(f"    pairs restricted to Cb-Cb > {FAR} A and |i-j| >= 6")
print(f"{'target':26s} {'|T|':>4s} {'npairs':>7s} {'loop coupling':>14s} {'null mean':>10s} "
      f"{'pct':>6s} {'p':>8s}")
W = len(LOOP)
windows = [list(range(s, s+W)) for s in range(1, L-W+2)]
rows = []
for tname, T in TARGETS.items():
    obs, npair = noncontact_coupling(LOOP, T)
    nulls = []
    for w in windows:
        if set(w) & set(LOOP): continue
        v, n = noncontact_coupling(w, T)
        if np.isfinite(v) and n >= 20: nulls.append(v)
    nulls = np.array(nulls)
    pct = 100*np.mean(nulls < obs); p = (np.sum(nulls >= obs)+1)/(len(nulls)+1)
    print(f"{tname:26s} {len(T):4d} {npair:7d} {obs:14.5f} {nulls.mean():10.5f} "
          f"{pct:6.1f} {p:8.4f}")
    rows.append(dict(target=tname, n_target=len(T), n_pairs=npair, loop=float(obs),
                     null_mean=float(nulls.mean()), percentile=float(pct), p=float(p),
                     n_null_windows=int(len(nulls))))

print(f"\n=== (b) do TUBA8's divergent positions coevolve as a module? ===")
def within(sel):
    vals=[]
    for a_ in range(len(sel)):
        for b_ in range(a_+1, len(sel)):
            i, j = sel[a_], sel[b_]
            if abs(i-j) < 6: continue
            vals.append(DI[i-1, j-1])
    return np.array(vals)
obs = within(DIV); 
rng = np.random.default_rng(3)
nulls = np.array([within(sorted(rng.choice(np.arange(1, L+1), size=len(DIV), replace=False))).mean()
                  for _ in range(20000)])
print(f"    {len(DIV)} divergent positions, {len(obs)} pairs with |i-j|>=6")
print(f"    mean within-set coupling = {obs.mean():.5f};  null {nulls.mean():.5f} "
      f"+/- {nulls.std():.5f}")
print(f"    p_enrich = {(np.sum(nulls>=obs.mean())+1)/20001:.4f}   "
      f"z = {(obs.mean()-nulls.mean())/nulls.std():+.2f}")
# restrict to non-contacting divergent pairs
obs_far=[]
for a_ in range(len(DIV)):
    for b_ in range(a_+1, len(DIV)):
        i,j = DIV[a_], DIV[b_]
        if abs(i-j) < 6: continue
        d = D[i-1,j-1]
        if np.isfinite(d) and d > FAR: obs_far.append(DI[i-1,j-1])
print(f"    restricted to spatially distant divergent pairs (>{FAR} A): n={len(obs_far)}, "
      f"mean {np.mean(obs_far):.5f}")

print(f"\n=== (c) strongest NON-CONTACTING couplings involving the loop or position 40 ===")
cands=[]
for i in list(LOOP)+[40]:
    for j in range(1, L+1):
        if abs(i-j) < 6: continue
        d = D[i-1, j-1]
        if not np.isfinite(d) or d <= FAR: continue
        cands.append((DI[i-1,j-1], i, j, d))
cands.sort(reverse=True)
print(f"    {'DI':>8s} {'i':>4s} {'j':>4s} {'Cb-Cb':>7s}  region of j")
REG = {d["canon_pos"]: d["region"] for d in div}
tb = set(TARGETS["TBCB interface (AF3)"])
for v,i,j,d in cands[:15]:
    tag = []
    if j in tb: tag.append("TBCB-iface")
    if j in ctx and ctx[j].get("d_nucleotide",99) < 5: tag.append("nucleotide")
    if j in ctx and ctx[j]["dSASA_beta"] > 5: tag.append("ab-iface")
    if j in REG: tag.append("TUBA8-divergent")
    print(f"    {v:8.4f} {i:4d} {j:4d} {d:7.1f}  {','.join(tag) if tag else '-'}")
json.dump(dict(loop_coupling=rows,
               module=dict(mean=float(obs.mean()), null_mean=float(nulls.mean()),
                           null_sd=float(nulls.std()), n_pairs=int(len(obs)),
                           p=float((np.sum(nulls>=obs.mean())+1)/20001)),
               far_pairs=[dict(DI=float(v), i=int(i), j=int(j), dist=float(d)) for v,i,j,d in cands[:40]]),
          open("results/dca/allostery.json","w"), indent=1)
