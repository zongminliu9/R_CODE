#!/usr/bin/env python3
"""AlphaFold pLDDT around the 32-45 loop, parsed directly from the B-factor column.

Built-in negative control: TUBA1A, TUBA1B, TUBA1C, TUBA3A and TUBA4A have an IDENTICAL 32-45
segment, so any spread among them measures model noise rather than sequence effect. TUBA8's
deviation is only interpretable if it exceeds that spread."""
import json, numpy as np
MODELS = {"TUBA8":"Q9JJZ2","TUBA1A":"P68369","TUBA1B":"P05213","TUBA4A":"P68368",
          "hTUBA8":"Q9NY65"}
LOOP = list(range(32, 46))

def plddt(path):
    per = {}
    for l in open(path):
        if not l.startswith("ATOM"): continue
        rid = int(l[22:26]); b = float(l[60:66])
        per.setdefault(rid, []).append(b)
    return {r: float(np.mean(v)) for r, v in per.items()}

print(f"{'isoform':9s} {'loop pLDDT':>11s} {'rest':>8s} {'delta':>7s} {'loop min':>9s} "
      f"{'windows below':>14s}")
out = {}
for iso, acc in MODELS.items():
    p = plddt(f"data/afdb/{acc}.pdb")
    loopv = [p[r] for r in LOOP if r in p]
    restv = [v for r, v in p.items() if r not in LOOP]
    keys = sorted(p)
    wins = np.array([np.mean([p[k] for k in keys[i:i+14]]) for i in range(len(keys)-13)])
    rank = int((wins < np.mean(loopv)).sum())
    print(f"{iso:9s} {np.mean(loopv):11.2f} {np.mean(restv):8.2f} "
          f"{np.mean(loopv)-np.mean(restv):+7.2f} {min(loopv):9.2f} "
          f"{f'{rank}/{len(wins)}':>14s}")
    out[iso] = dict(loop=float(np.mean(loopv)), rest=float(np.mean(restv)),
                    delta=float(np.mean(loopv)-np.mean(restv)), loop_min=float(min(loopv)),
                    rank=rank, n_windows=len(wins),
                    per_residue={str(r): round(p[r],2) for r in LOOP if r in p})

ctrl = [out[i]["loop"] for i in ("TUBA1A","TUBA1B","TUBA4A")]
print(f"\ncontrol spread (TUBA1A/1B/4A, identical 32-45 sequence): "
      f"{min(ctrl):.2f}-{max(ctrl):.2f}  (range {max(ctrl)-min(ctrl):.2f})")
print(f"TUBA8 mouse loop pLDDT {out['TUBA8']['loop']:.2f}, human TUBA8 {out['hTUBA8']['loop']:.2f}")
d = out["TUBA8"]["loop"] - np.mean(ctrl)
print(f"TUBA8 deviation from control mean: {d:+.2f}  -> "
      f"{'EXCEEDS' if abs(d) > (max(ctrl)-min(ctrl)) else 'WITHIN'} the control spread")

print("\nper-residue pLDDT across the loop (mouse):")
print("  pos  " + "".join(f"{r:>7d}" for r in LOOP))
for iso in ("TUBA1A","TUBA1B","TUBA4A","TUBA8"):
    pr = out[iso]["per_residue"]
    print(f"  {iso:7s}" + "".join(f"{pr.get(str(r),float('nan')):7.1f}" for r in LOOP))
json.dump(out, open("results/structure/plddt_loop32_45.json","w"), indent=1)
