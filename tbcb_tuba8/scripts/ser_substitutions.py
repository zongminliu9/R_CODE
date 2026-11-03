#!/usr/bin/env python3
"""The scattered surface substitutions, and the apparent excess of serines.

Questions: is the shift toward Ser/Thr real or an impression? where are those positions, and
are they clustered? do they change the surface electrostatics? and -- using the DCA couplings --
do they coevolve with each other specifically, over and above the whole divergent set?
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import binomtest, mannwhitneyu
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx

div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
ctx = {r["canon_pos"]: r for r in json.load(open("results/structure/alpha_residue_structural_context_9M1M.json"))}
DI  = np.load("results/dca/dca_APC_full.npy").astype(np.float64); L = DI.shape[0]
idx = json.load(open("data/meta/sequence_index.json"))
a1a = idx["mouse_Tuba1a"]["sequence"]
rng = np.random.default_rng(77)

print("=== substitution spectrum TUBA1A -> TUBA8 (45 positions) ===")
to_counts = {}; from_counts = {}; pairs = []
for d in div:
    f, t = d["TUBA1A"], d["TUBA8"]
    if t == "-" or f == "-": continue
    to_counts[t] = to_counts.get(t, 0)+1
    from_counts[f] = from_counts.get(f, 0)+1
    pairs.append((d["canon_pos"], f, t, d["region"], d["tuba8_unique"]))
print("  gained (TUBA8 side): " + "  ".join(f"{k}:{v}" for k,v in sorted(to_counts.items(), key=lambda x:-x[1])))
print("  lost   (TUBA1A side):" + "  ".join(f"{k}:{v}" for k,v in sorted(from_counts.items(), key=lambda x:-x[1])))

n_sub = len(pairs)
ser_thr_gained = [p for p in pairs if p[2] in "ST"]
ser_gained     = [p for p in pairs if p[2] == "S"]
print(f"\n  substitutions to Ser or Thr: {len(ser_thr_gained)}/{n_sub} "
      f"({100*len(ser_thr_gained)/n_sub:.1f}%)   to Ser alone: {len(ser_gained)}")
print("   " + "  ".join(f"{p[1]}{p[0]}{p[2]}" for p in ser_thr_gained))

# expectation: background frequency of S+T among all residues of TUBA1A
bg_st = (a1a.count("S")+a1a.count("T"))/len(a1a)
bg_s  = a1a.count("S")/len(a1a)
bt = binomtest(len(ser_thr_gained), n_sub, bg_st, alternative="greater")
bs = binomtest(len(ser_gained), n_sub, bg_s, alternative="greater")
print(f"  background S+T frequency in TUBA1A = {100*bg_st:.1f}%  -> expected "
      f"{n_sub*bg_st:.1f}, observed {len(ser_thr_gained)},  binomial p = {bt.pvalue:.4f}")
print(f"  background S   frequency = {100*bg_s:.1f}%  -> expected {n_sub*bg_s:.1f}, "
      f"observed {len(ser_gained)},  p = {bs.pvalue:.4f}")
# net Ser/Thr change (gained minus lost)
lost_st = sum(1 for p in pairs if p[1] in "ST")
print(f"  NET change: TUBA8 gains {len(ser_thr_gained)} S/T and loses {lost_st} -> "
      f"net {len(ser_thr_gained)-lost_st:+d} hydroxyl side chains at divergent positions")

print("\n=== charge and polarity shift ===")
chg = lambda s: sum(1 for c in s if c in "KR") - sum(1 for c in s if c in "DE")
a8 = idx["mouse_Tuba8"]["sequence"]
print(f"  whole protein net charge: TUBA1A {chg(a1a):+d}  TUBA8 {chg(a8):+d}")
dq = sum((1 if t in "KR" else -1 if t in "DE" else 0) - (1 if f in "KR" else -1 if f in "DE" else 0)
         for _, f, t, _, _ in pairs)
print(f"  charge change contributed by the 45 divergent positions: {dq:+d}")
arom = lambda s: sum(1 for c in s if c in "FWY")
print(f"  aromatics: TUBA1A {arom(a1a)}  TUBA8 {arom(a8)}   "
      f"(divergent-position aromatics gained: "
      f"{sum(1 for _,f,t,_,_ in pairs if t in 'FWY' and f not in 'FWY')}, "
      f"lost: {sum(1 for _,f,t,_,_ in pairs if f in 'FWY' and t not in 'FWY')})")

print("\n=== where are the gained Ser/Thr, and are they buried or exposed? ===")
print(f"  {'pos':>5s} {'sub':>6s} {'relSASA':>8s} {'dSASA_cof':>10s} {'region'}")
for p, f, t, reg, uq in ser_thr_gained:
    c = ctx.get(p)
    rs = f"{c['relSASA_mono']:.3f}" if c else "n/a"
    dc = f"{c['dSASA_cofactor']:.1f}" if c else "n/a"
    print(f"  {p:5d} {f+str(p)+t:>6s} {rs:>8s} {dc:>10s} {reg}")

print("\n=== spatial clustering of gained Ser/Thr (9M1M chain a) ===")
arr = pdbx.get_structure(pdbx.CIFFile.read("data/pdb/9M1M.cif"), model=1, use_author_fields=True)
A = arr[(arr.chain_id=="a") & struc.filter_amino_acids(arr)]
ca = {int(r): A[(A.res_id==r)&(A.atom_name=="CA")].coord[0] for r in np.unique(A.res_id)
      if A[(A.res_id==r)&(A.atom_name=="CA")].array_length()}
st_pos = [p for p,_,_,_,_ in ser_thr_gained if p in ca]
if len(st_pos) >= 3:
    def mean_pair_dist(sel):
        d=[]
        for i in range(len(sel)):
            for j in range(i+1,len(sel)):
                d.append(np.linalg.norm(ca[sel[i]]-ca[sel[j]]))
        return np.mean(d)
    obs = mean_pair_dist(st_pos)
    keys = sorted(ca)
    null = np.array([mean_pair_dist(list(rng.choice(keys, size=len(st_pos), replace=False)))
                     for _ in range(20000)])
    print(f"  {len(st_pos)} gained S/T resolved; mean pairwise CA distance {obs:.2f} A")
    print(f"  random null {null.mean():.2f} +/- {null.std():.2f} A; "
          f"p(clustered) = {(np.sum(null<=obs)+1)/20001:.4f}")

print("\n=== do the gained Ser/Thr coevolve with each other specifically? ===")
def within(sel):
    v=[]
    for i in range(len(sel)):
        for j in range(i+1,len(sel)):
            if abs(sel[i]-sel[j])<6: continue
            v.append(DI[sel[i]-1, sel[j]-1])
    return np.array(v)
sel=[p for p,_,_,_,_ in ser_thr_gained]
o=within(sel)
if len(o):
    allsel=[p for p,_,_,_,_ in pairs]
    nulls=np.array([within(list(rng.choice(allsel, size=len(sel), replace=False))).mean()
                    for _ in range(20000)])
    print(f"  mean coupling among gained S/T = {o.mean():.5f} ({len(o)} pairs)")
    print(f"  null = random subsets of the SAME divergent set: {nulls.mean():.5f} "
          f"+/- {nulls.std():.5f}; z = {(o.mean()-nulls.mean())/nulls.std():+.2f}; "
          f"p = {(np.sum(nulls>=o.mean())+1)/20001:.4f}")
