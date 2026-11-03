#!/usr/bin/env python3
"""Is TUBA8 uniquely destabilised at BURIED positions?

The kinetic model's surviving requirement is a shorter lifetime for the chaperone-free
quasi-native state (k_degF), which is a FOLDED-STATE STABILITY property, not an interface
property. Surface substitutions cost little folding free energy; buried ones cost a lot. So the
question is not how many substitutions an isoform carries but how much buried packing they
disturb.

Proxy (and it is a proxy, not a free energy): for each substitution,
    cost_i = burial_i * [ dG_transfer(from) - dG_transfer(to) ]_+        hydrophobic burial lost
           + burial_i * volume_penalty(|dV|)                            packing strain
           + burial_i * polar_burial(to)                                buried polar desolvation
with burial_i = 1 - relSASA_i taken from the experimental 9M1M coordinates. Every isoform is
scored by the SAME proxy against the SAME reference, so the comparison between isoforms is
meaningful even though an absolute value is not.
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_isoforms import nw_affine

# Fauchere-Pliska octanol-water transfer free energies (kcal/mol, positive = hydrophobic)
FP = {'W':2.25,'I':1.80,'F':1.79,'L':1.70,'C':1.54,'M':1.23,'V':1.22,'Y':0.96,'P':0.72,
      'A':0.31,'T':0.26,'H':0.13,'G':0.00,'S':-0.04,'Q':-0.22,'N':-0.60,'E':-0.64,
      'D':-0.77,'K':-0.99,'R':-1.01}
# residue volumes (A^3, Zamyatnin)
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,
       'H':153.2,'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,
       'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
POLAR = set("STNQDEKRH")

idx = json.load(open("data/meta/sequence_index.json"))
ctx = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
ref = idx["mouse_Tuba1a"]["sequence"]
ISO = ["Tuba8","Tuba1b","Tuba1c","Tuba3a","Tuba4a"]

def subs_vs_ref(q):
    A, B = nw_affine(ref, q); i = j = 0; out = []
    for x, y in zip(A, B):
        if x != "-": i += 1
        if y != "-": j += 1
        if x != "-" and y != "-" and x != y: out.append((i, x, y))
    return out

def burial(p):
    c = ctx.get(p)
    if c is None: return None
    return max(0.0, min(1.0, 1.0 - c["relSASA_mono"]))

def cost(f, t, b):
    dh = max(0.0, FP.get(f,0.0) - FP.get(t,0.0))          # hydrophobicity lost
    dv = abs(VOL.get(f,120.0) - VOL.get(t,120.0))/100.0   # packing mismatch, per 100 A^3
    pp = 1.0 if (t in POLAR and f not in POLAR) else 0.0  # newly buried polar group
    return b*dh + 0.6*b*dv + 0.8*b*pp

print(f"{'isoform':9s} {'n_sub':>6s} {'n_resolved':>11s} {'n_buried(b>0.8)':>16s} "
      f"{'sum cost':>9s} {'mean cost':>10s} {'max single':>11s}")
tab = {}
for g in ISO:
    q = idx[f"mouse_{g}"]["sequence"]
    S = subs_vs_ref(q)
    rows = []
    for p, f, t in S:
        b = burial(p)
        if b is None: continue
        rows.append((p, f, t, b, cost(f, t, b)))
    nb = sum(1 for r in rows if r[3] > 0.8)
    tot = sum(r[4] for r in rows)
    print(f"{g.upper():9s} {len(S):6d} {len(rows):11d} {nb:16d} {tot:9.3f} "
          f"{tot/max(len(rows),1):10.3f} {max((r[4] for r in rows), default=0):11.3f}")
    tab[g] = rows

print("\n=== TUBA8: substitutions ranked by the burial-weighted cost proxy ===")
print(f"  {'pos':>5s} {'sub':>7s} {'burial':>7s} {'relSASA':>8s} {'dFP':>6s} {'dVol':>7s} "
      f"{'cost':>6s}  region")
div = {d["canon_pos"]: d["region"] for d in json.load(open("results/alignment/tuba8_divergent_positions.json"))}
for p, f, t, b, c in sorted(tab["Tuba8"], key=lambda r: -r[4])[:14]:
    print(f"  {p:5d} {f+str(p)+t:>7s} {b:7.3f} {1-b:8.3f} "
          f"{FP.get(f,0)-FP.get(t,0):6.2f} {VOL.get(f,0)-VOL.get(t,0):7.1f} {c:6.3f}  {div.get(p,'')}")

print("\n=== is TUBA8's buried-substitution load an outlier? ===")
tot8 = sum(r[4] for r in tab["Tuba8"])
others = {g: sum(r[4] for r in tab[g]) for g in ISO if g != "Tuba8"}
print(f"  TUBA8 total cost {tot8:.3f}")
for g, v in sorted(others.items(), key=lambda x:-x[1]):
    print(f"  {g.upper():8s} {v:8.3f}   TUBA8/{g.upper()} = {tot8/max(v,1e-9):5.2f}x")
# per-substitution normalised: does TUBA8 just have more substitutions, or worse ones?
print("\n  normalised by substitution count (mean cost per substitution):")
for g in ISO:
    rows = tab[g]
    print(f"    {g.upper():8s} {sum(r[4] for r in rows)/max(len(rows),1):6.3f} "
          f"over {len(rows)} resolved substitutions")
# permutation: given TUBA8's 39 resolved positions, is the cost higher than random substitutions
# placed at the same positions with the same amino-acid composition?
rng = np.random.default_rng(42)
obs = tot8
comp = [(r[1], r[2]) for r in tab["Tuba8"]]
positions = [r[0] for r in tab["Tuba8"]]
allpos = [p for p in ctx]
null = []
for _ in range(20000):
    ps = rng.choice(allpos, size=len(comp), replace=False)
    null.append(sum(cost(f, t, burial(p)) for p, (f, t) in zip(ps, comp)))
null = np.array(null)
print(f"\n  keeping TUBA8's substitution identities but RANDOMISING their positions:")
print(f"    observed {obs:.3f}; null {null.mean():.3f} +/- {null.std():.3f}; "
      f"z = {(obs-null.mean())/null.std():+.2f}; "
      f"p(higher) = {(np.sum(null>=obs)+1)/20001:.4f}")
json.dump({g: [[int(p),f,t,float(b),float(c)] for p,f,t,b,c in tab[g]] for g in ISO},
          open("results/structure/burial_cost.json","w"), indent=1)
