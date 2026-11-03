#!/usr/bin/env python3
"""Map the 45 TUBA8-divergent canonical positions onto experimental structures and ask, without
bias toward pre-named residues, WHERE they sit: buried core / nucleotide pocket / alpha-beta
intradimer interface / cofactor (TBCD,TBCE,TBCC,Arl2) contact surface / free surface / C-term tail.

Structures: 1JFF (canonical ab-tubulin numbering) and 9M1M (TBC-DEC-Arl2-ab-tubulin, 2.21 A).
Chains are identified automatically by sequence identity to the UniProt references, so no
hard-coded chain letters can silently go stale.
"""
import json, csv, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
# Tien et al. 2013 theoretical max SASA (A^2), Gly-X-Gly, for relative-SASA normalisation
MAXSASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
           'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}

seqidx = json.load(open("data/meta/sequence_index.json"))
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
DIVPOS = {d["canon_pos"]: d for d in div}

def load(cif):
    f = pdbx.CIFFile.read(cif)
    arr = pdbx.get_structure(f, model=1, use_author_fields=True, include_bonds=False)
    return arr

def chain_seq(arr, ch):
    sub = arr[(arr.chain_id == ch) & struc.filter_amino_acids(arr)]
    if sub.array_length() == 0: return "", []
    ids, names = [], []
    seen = set()
    for rid, rn in zip(sub.res_id, sub.res_name):
        if (rid,) in seen: continue
        seen.add((rid,)); ids.append(int(rid)); names.append(rn)
    return "".join(AA3.get(n, "X") for n in names), ids

def best_match(seq, refs):
    """crude but decisive: fraction of reference 12-mers present in the chain sequence."""
    best, bestsc = None, 0.0
    for name, rs in refs.items():
        kmers = {rs[i:i+12] for i in range(0, len(rs)-12, 6)}
        if not kmers: continue
        sc = sum(1 for k in kmers if k in seq)/len(kmers)
        if sc > bestsc: best, bestsc = name, sc
    return best, bestsc

REFS = {k: v["sequence"] for k, v in seqidx.items() if k.startswith(("mouse_","human_"))}

def analyse(cif, tag):
    arr = load(cif)
    prot = arr[struc.filter_amino_acids(arr)]
    chains = [c for c in np.unique(prot.chain_id)]
    ident = {}
    for ch in chains:
        s, ids = chain_seq(prot, ch)
        if len(s) < 40: continue
        nm, sc = best_match(s, REFS)
        ident[ch] = (nm, round(sc,2), len(s), ids[0] if ids else None, ids[-1] if ids else None)
    print(f"\n=== {tag} ({os.path.basename(cif)}) chain identification ===")
    for ch,(nm,sc,n,lo,hi) in sorted(ident.items()):
        print(f"  chain {ch:>2s}: {nm or '?':20s} match={sc:4.2f}  {n:4d} res  auth {lo}-{hi}")
    return arr, prot, ident

if __name__ == "__main__":
    for cif, tag in [("data/pdb/1JFF.cif","ab-tubulin reference"),
                     ("data/pdb/9M1M.cif","TBC-DEC-Arl2-ab-tubulin")]:
        analyse(cif, tag)
