#!/usr/bin/env python3
"""Read AlphaFold 3 confidence metrics for the uploaded predictions.

The decisive number for a two-chain prediction is ipTM: it scores the INTERFACE, not the folds.
Below ~0.6 the interface placement is unreliable and nothing downstream should be built on it;
above ~0.8 it is usually right. pTM scores the whole complex and is inflated by the two
well-folded chains, so it must not be used in ipTM's place.
"""
import json, glob, os, sys
B = "/home/user/uploads_fold"
JOBS = {
 "TBCB+TUBA8":  f"{B}/e50dbe95-fold_tbcb_tuba8",
 "TBCB+TUBA1A": f"{B}/88398339-fold_tbcb_tuba1a",
 "TUBA1A":      f"{B}/a98a190e-fold_tuba1a",
 "TUBA1B":      f"{B}/96c4e58d-fold_tuba1b",
 "TUBA1A_31_50":f"{B}/ce252fa4-fold_tuba1a_31_50",
}
out = {}
for name, d in JOBS.items():
    print(f"\n=== {name} ===")
    # job request tells us what was actually submitted
    jr = glob.glob(f"{d}/*job_request.json")
    if jr:
        j = json.load(open(jr[0]))
        seqs = j[0].get("sequences", j[0]) if isinstance(j, list) else j.get("sequences", [])
        for s in (seqs if isinstance(seqs, list) else []):
            pc = s.get("proteinChain") or s.get("protein") or {}
            sq = pc.get("sequence", "")
            if sq: print(f"  submitted chain: {len(sq)} aa  {sq[:40]}...{sq[-12:]}")
        if j and isinstance(j, list) and "modelSeeds" in j[0]:
            print(f"  seeds: {j[0]['modelSeeds']}")
    rows = []
    for f in sorted(glob.glob(f"{d}/*summary_confidences_*.json")):
        c = json.load(open(f))
        k = os.path.basename(f).split("_")[-1].split(".")[0]
        rows.append((k, c))
    if not rows:
        print("  no summary_confidences files"); continue
    keys = ["ptm","iptm","ranking_score","fraction_disordered","has_clash",
            "num_recycles"]
    hdr = f"  {'model':>6s}" + "".join(f"{k:>20s}" for k in keys if k in rows[0][1])
    print(hdr)
    for k, c in rows:
        line = f"  {k:>6s}"
        for kk in keys:
            if kk in c:
                v = c[kk]
                line += f"{v:>20.4f}" if isinstance(v,(int,float)) else f"{str(v):>20s}"
        print(line)
    # chain-pair iptm matrix (the per-interface score)
    c0 = rows[0][1]
    for field in ("chain_pair_iptm","chain_pair_pae_min","chain_iptm","chain_ptm"):
        if field in c0:
            print(f"  {field}:")
            v = c0[field]
            if isinstance(v, list) and v and isinstance(v[0], list):
                for i, r in enumerate(v):
                    print(f"     {i}: " + " ".join(f"{x:7.3f}" if isinstance(x,(int,float)) else str(x) for x in r))
            else:
                print(f"     {v}")
    out[name] = {k: c for k, c in rows}
json.dump(out, open("results/af3/confidences.json","w"), indent=1)
