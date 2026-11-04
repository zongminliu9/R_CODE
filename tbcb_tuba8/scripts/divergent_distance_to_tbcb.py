#!/usr/bin/env python3
"""How far is each TUBA8 divergent position from TBCB in the AF3 complexes?

This converts the intermediate-state objection into a number. If a divergent position is 40 A
from TBCB in the native-state model, no plausible partial-folding rearrangement brings it into
contact and it is irrelevant to binding in any state. If it is 8-15 A away, a modest
repositioning of an unpacked domain could. The distance distribution therefore says which of
TUBA8's 38 core substitutions could matter for TBCB capture of a folding intermediate, and which
cannot matter in any state.

Also checks the AF3 model against an independent criterion: which TBCB-interface residues of
alpha-tubulin are the ones beta-tubulin uses in the mature heterodimer. A predicted interface
that sits on another partner's binding site is the signature of a misdocked model.
"""
import glob, json, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx

B = "/home/user/uploads_fold"
JOBS = {"TUBA8": (f"{B}/e50dbe95-fold_tbcb_tuba8", 449),
        "TUBA1A": (f"{B}/88398339-fold_tbcb_tuba1a", 451)}
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
ctx = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
af3 = json.load(open("results/af3/interface_aligned.json"))
DOM = lambda p: ("N-GTPase" if p<=205 else "intermediate" if p<=381 else
                 "C-term dom" if p<=440 else "tail")

# distance of every alpha residue to TBCB, averaged over the 5 models, in the TUBA8 complex
d, L = JOBS["TUBA8"]
cmap = {int(k): v for k, v in af3["TUBA8"]["cmap"].items()}   # raw -> canonical
inv  = {v: k for k, v in cmap.items()}
per = {}
for m in sorted(glob.glob(f"{d}/*_model_*.cif")):
    arr = pdbx.get_structure(pdbx.CIFFile.read(m), model=1, use_author_fields=True)
    arr = arr[struc.filter_amino_acids(arr)]
    chs = {c: int((arr.chain_id==c).sum()) for c in np.unique(arr.chain_id)}
    chT = min(chs, key=chs.get); chA = max(chs, key=chs.get)
    A = arr[arr.chain_id==chA]; T = arr[arr.chain_id==chT]
    dm = np.linalg.norm(A.coord[:,None,:]-T.coord[None,:,:], axis=-1)
    mn = dm.min(axis=1)
    for rid in np.unique(A.res_id):
        v = float(mn[A.res_id==rid].min())
        per.setdefault(int(rid), []).append(v)
dist = {r: float(np.mean(v)) for r, v in per.items()}
sd   = {r: float(np.std(v))  for r, v in per.items()}

print("=== distance from TBCB (min heavy-atom, mean over 5 AF3 models) ===")
print(f"{'canon':>6s} {'sub':>8s} {'domain':>13s} {'dist(A)':>9s} {'sd':>5s} "
      f"{'relSASA':>8s}  verdict")
rows = []
for dd in sorted(div, key=lambda x: x["canon_pos"]):
    p = dd["canon_pos"]; raw = inv.get(p)
    if raw is None or raw not in dist: continue
    v = dist[raw]
    verdict = ("CONTACT" if v < 4.5 else
               "near (<10A): reachable by modest rearrangement" if v < 10 else
               "intermediate (10-20A)" if v < 20 else
               "far (>20A): cannot contact TBCB in any folded-domain arrangement")
    rs = ctx[p]["relSASA_mono"] if p in ctx else float("nan")
    print(f"{p:6d} {dd['TUBA1A']+str(p)+dd['TUBA8']:>8s} {DOM(p):>13s} {v:9.1f} "
          f"{sd[raw]:5.1f} {rs:8.3f}  {verdict}")
    rows.append(dict(canon=p, sub=f"{dd['TUBA1A']}{p}{dd['TUBA8']}", domain=DOM(p),
                     dist=v, sd=sd[raw], relSASA=None if np.isnan(rs) else rs))
arr_d = np.array([r["dist"] for r in rows])
print(f"\n  {len(rows)} divergent positions mapped")
for lo, hi, lab in ((0,4.5,"in contact"),(4.5,10,"within 10 A"),(10,20,"10-20 A"),
                    (20,1e9,"beyond 20 A")):
    n = int(((arr_d>=lo)&(arr_d<hi)).sum())
    print(f"    {lab:14s} {n:3d}  ({100*n/len(rows):4.1f}%)")
core = [r for r in rows if r["domain"] != "tail"]
print(f"\n  core (non-tail) divergent positions: {len(core)}")
print(f"    median distance to TBCB {np.median([r['dist'] for r in core]):.1f} A, "
      f"min {min(r['dist'] for r in core):.1f} A "
      f"({[r['sub'] for r in core if r['dist']==min(x['dist'] for x in core)]})")
within10 = [r for r in core if r["dist"] < 10]
print(f"    core positions within 10 A of TBCB: "
      f"{[r['sub'] for r in within10] if within10 else 'NONE'}")

print("\n=== does the AF3 TBCB interface sit on beta-tubulin's binding site? ===")
tbcb = set(af3["TUBA8"]["robust_canon"]) | set(af3["TUBA1A"]["robust_canon"])
beta = {p for p,r in ctx.items() if r["dSASA_beta"] > 5}
ov = sorted(tbcb & beta)
print(f"  beta-tubulin interface (9M1M): {sorted(beta)}")
print(f"  shared with the AF3 TBCB interface: {ov}  ({len(ov)}/{len(beta)} of beta's site)")
if ov:
    dS = {p: ctx[p]["dSASA_beta"] for p in ov}
    print(f"  those residues bury {sum(dS.values()):.0f} A^2 against beta-tubulin, "
          f"i.e. {100*sum(dS.values())/sum(ctx[p]['dSASA_beta'] for p in beta):.0f}% "
          f"of the intradimer interface")
json.dump(dict(distances=rows, beta_overlap=ov), 
          open("results/af3/divergent_distance_to_tbcb.json","w"), indent=1)
