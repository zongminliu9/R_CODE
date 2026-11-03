#!/usr/bin/env python3
"""Compare the two predicted TBCB interfaces in ALIGNED coordinates.

TUBA8 is 449 aa and TUBA1A 451 aa, so raw residue numbers are offset by 2 in the C-terminal
tail. Comparing them directly, as the first pass did, manufactures fake differences there. Here
every interface residue is mapped to canonical TUBA1A numbering through the global alignment
before the two are compared.

Also: which part of TBCB engages the tail? TBCB has an N-terminal ubiquitin-like domain and a
C-terminal CAP-Gly domain at Q9D1E6 residues 134-231. If the tail is read by the CAP-Gly domain
that is the canonical mechanism; if not, it is something new.
"""
import json, glob, os, warnings, sys
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_isoforms import nw_affine

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
B = "/home/user/uploads_fold"
JOBS = {"TUBA8": f"{B}/e50dbe95-fold_tbcb_tuba8", "TUBA1A": f"{B}/88398339-fold_tbcb_tuba1a"}
CAPGLY = set(range(134, 232)); UBL = set(range(1, 100))
idx = json.load(open("data/meta/sequence_index.json"))
ref = idx["mouse_Tuba1a"]["sequence"]
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
DIV = {d["canon_pos"] for d in div}

def canon_map(qseq):
    """query residue index (1-based) -> canonical TUBA1A position"""
    A, Bq = nw_affine(ref, qseq); i = j = 0; m = {}
    for x, y in zip(A, Bq):
        if x != "-": i += 1
        if y != "-": j += 1
        if x != "-" and y != "-": m[j] = i
    return m

def analyse(d, iso):
    models = sorted(glob.glob(f"{d}/*_model_*.cif"))
    acc = {}; tb_acc = {}; pair_acc = {}
    qseq = None
    for m in models:
        arr = pdbx.get_structure(pdbx.CIFFile.read(m), model=1, use_author_fields=True)
        arr = arr[struc.filter_amino_acids(arr)]
        chs = {c: (arr[arr.chain_id==c][arr[arr.chain_id==c].atom_name=="CA"]).array_length()
               for c in np.unique(arr.chain_id)}
        chT = min(chs, key=chs.get); chA = max(chs, key=chs.get)
        A = arr[arr.chain_id == chA]; T = arr[arr.chain_id == chT]
        if qseq is None:
            ca = A[A.atom_name == "CA"]; o = np.argsort(ca.res_id)
            qseq = "".join(AA3.get(n,"X") for n in ca.res_name[o])
        dm = np.linalg.norm(A.coord[:,None,:]-T.coord[None,:,:], axis=-1)
        # alpha residues in contact
        for i in np.where(dm.min(axis=1) < 4.5)[0]:
            acc[int(A.res_id[i])] = acc.get(int(A.res_id[i]), 0)+1
        # TBCB residues in contact
        for j in np.where(dm.min(axis=0) < 4.5)[0]:
            tb_acc[int(T.res_id[j])] = tb_acc.get(int(T.res_id[j]), 0)+1
        # residue-residue pairs
        ii, jj = np.where(dm < 4.5)
        for i, j in zip(ii, jj):
            k = (int(A.res_id[i]), int(T.res_id[j]))
            pair_acc[k] = pair_acc.get(k, 0)+1
    nmod = len(models)
    return qseq, {r:c for r,c in acc.items()}, tb_acc, pair_acc, nmod

OUT = {}
for iso, d in JOBS.items():
    qseq, acc, tb, pairs, nmod = analyse(d, iso)
    cm = canon_map(qseq)
    robust = sorted(r for r,c in acc.items() if c >= 4)
    canon = sorted({cm[r] for r in robust if r in cm})
    OUT[iso] = dict(qlen=len(qseq), robust_raw=robust, robust_canon=canon,
                    tbcb_contacts=sorted(r for r,c in tb.items() if c >= 4),
                    pairs={f"{a}|{b}": c for (a,b),c in pairs.items() if c >= 4},
                    cmap={str(k): v for k,v in cm.items()})
    print(f"\n=== TBCB + {iso} ===")
    print(f"  interface in CANONICAL TUBA1A numbering ({len(canon)} residues):")
    print("    " + " ".join(f"{p}{'*' if p in DIV else ''}" for p in canon))
    tbc = OUT[iso]["tbcb_contacts"]
    ncap = [r for r in tbc if r in CAPGLY]; nubl = [r for r in tbc if r in UBL]
    nmid = [r for r in tbc if r not in CAPGLY and r not in UBL]
    print(f"  TBCB side: {len(tbc)} residues -> CAP-Gly(134-231): {len(ncap)}, "
          f"UBL(1-99): {len(nubl)}, linker/other: {len(nmid)}")
    print(f"    CAP-Gly contacts: {ncap}")
    print(f"    UBL contacts:     {nubl}")
    print(f"    linker/other:     {nmid}")

c8 = set(OUT["TUBA8"]["robust_canon"]); c1 = set(OUT["TUBA1A"]["robust_canon"])
print(f"\n=== aligned comparison ===")
print(f"  TUBA8 interface {len(c8)} residues; TUBA1A {len(c1)}; shared {len(c8&c1)}")
print(f"  TUBA8-only (canonical): {sorted(c8-c1)}")
print(f"  TUBA1A-only (canonical): {sorted(c1-c8)}")
print(f"  divergent positions at interface: TUBA8 {sorted(c8&DIV)}  TUBA1A {sorted(c1&DIV)}")

# permutation test: are divergent positions enriched at the TBCB interface?
rng = np.random.default_rng(5)
allpos = np.arange(1, 452)
obs = len(c8 & DIV); k = len(DIV)
mask = np.isin(allpos, sorted(c8))
null = np.array([mask[rng.choice(len(allpos), size=k, replace=False)].sum() for _ in range(200000)])
print(f"\n  permutation test (TBCB interface, canonical 1-451, 2e5 draws):")
print(f"    observed divergent-at-interface = {obs}; expected {null.mean():.2f}; "
      f"p_enrich = {(np.sum(null>=obs)+1)/200001:.4f}; p_deplete = {(np.sum(null<=obs)+1)/200001:.4f}")
# split: tail vs core
TAIL = set(range(437, 452))
dt = [p for p in DIV if p in TAIL]; dc = [p for p in DIV if p not in TAIL]
print(f"    divergent positions in the tail (437-451): {sorted(dt)}; "
      f"at interface: {sorted(set(dt)&c8)}  ({len(set(dt)&c8)}/{len(dt)})")
print(f"    divergent positions in the core (<437): {len(dc)}; "
      f"at interface: {sorted(set(dc)&c8)}  ({len(set(dc)&c8)}/{len(dc)})")
json.dump(OUT, open("results/af3/interface_aligned.json","w"), indent=1)
