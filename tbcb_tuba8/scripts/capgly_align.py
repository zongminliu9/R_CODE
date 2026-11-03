#!/usr/bin/env python3
"""Align the CAP-Gly domains of mouse TBCB (1WHG), p150Glued (2HL3, which carries the bound
EEY tripeptide) and T. brucei TBCB (4B6M, 1.59 A), so the experimental EEY binding mode can be
transplanted into the mouse TBCB groove.

The CAP-Gly fold is recognised by its GKNDG motif, which forms one wall of the groove that
reads the -EEY/F-COO- terminus. We locate that motif in each domain explicitly rather than
trusting a sequence alignment blindly."""
import json, sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_isoforms import nw_affine

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}

def chain_ca(cif, chain):
    f = pdbx.CIFFile.read(cif)
    arr = pdbx.get_structure(f, model=1, use_author_fields=True)
    sel = arr[(arr.chain_id == chain) & struc.filter_amino_acids(arr)]
    ca  = sel[sel.atom_name == "CA"]
    o   = np.argsort(ca.res_id)
    seq = "".join(AA3.get(n, "X") for n in ca.res_name[o])
    return ca.coord[o].astype(float), ca.res_id[o].astype(int), seq, arr

import re
TARGETS = [("2HL3","A","p150Glued CAP-Gly (holds EEY peptide)"),
           ("1WHG","A","mouse TBCB CAP-Gly (NMR)"),
           ("4B6M","A","T. brucei TBCB CAP-Gly 1.59 A")]
info = {}
for pid, ch, desc in TARGETS:
    coord, rid, seq, arr = chain_ca(f"data/pdb/{pid}.cif", ch)
    # locate the GKNDG-type motif: G-x-x-D-G with a preceding hydrophobic/aromatic stretch
    hits = [m.start() for m in re.finditer(r"G[KR].[DN]G", seq)]
    hits2 = [m.start() for m in re.finditer(r"GKNDG|GKGDG|GKNDA", seq)]
    info[pid] = dict(coord=coord, rid=rid, seq=seq, arr=arr, chain=ch, desc=desc)
    print(f"{pid} ch{ch}  {len(seq):3d} res  auth {rid[0]}-{rid[-1]}  {desc}")
    print(f"   seq: {seq}")
    for h in sorted(set(hits+hits2)):
        print(f"   motif '{seq[h:h+5]}' at seq idx {h} -> auth {rid[h]}-{rid[h+4]}")
json.dump({k: dict(seq=v["seq"], rid=v["rid"].tolist(), chain=v["chain"], desc=v["desc"])
           for k, v in info.items()}, open("results/structure/capgly_domains.json","w"), indent=1)

# ---- pairwise sequence alignment to define equivalent positions
print("\n=== pairwise CAP-Gly sequence alignment ===")
def aln_pairs(sA, ridA, sB, ridB):
    A, B = nw_affine(sA, sB)
    i = j = 0; pairs = []
    for x, y in zip(A, B):
        if x != "-": i += 1
        if y != "-": j += 1
        if x != "-" and y != "-": pairs.append((ridA[i-1], ridB[j-1], x, y))
    ident = sum(1 for _,_,x,y in pairs if x == y)
    return pairs, ident
for a, b in [("2HL3","1WHG"), ("2HL3","4B6M"), ("1WHG","4B6M")]:
    pairs, ident = aln_pairs(info[a]["seq"], info[a]["rid"], info[b]["seq"], info[b]["rid"])
    print(f"  {a} vs {b}: {len(pairs)} aligned, {ident} identical ({100*ident/len(pairs):.1f}%)")
    json.dump([[int(p[0]),int(p[1]),p[2],p[3]] for p in pairs],
              open(f"results/structure/aln_{a}_{b}.json","w"))
