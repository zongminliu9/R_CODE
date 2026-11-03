#!/usr/bin/env python3
"""The kinetic model requires a difference between TBCB capture (kon_B) and direct TBCE capture
(kon_Edir). TBCE also carries a CAP-Gly domain, so both cofactors read the same -EEY/F tail.
That makes differential reading of the TUBA8 tail by TBCB-CAP-Gly versus TBCE-CAP-Gly a
structurally concrete candidate for exactly the asymmetry the model needs.

Locate TBCE's CAP-Gly domain in sequence and in the 9M1M coordinates, where TBCE chain E is
resolved over 1-527."""
import json, re, sys, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_isoforms import nw_affine

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
idx = json.load(open("data/meta/sequence_index.json"))

print("=== GKNDG-type motif search in the cofactor sequences ===")
for key in ("mouse_Tbcb","human_TBCB","mouse_Tbce","human_TBCE"):
    s = idx[key]["sequence"]
    hits = sorted({m.start() for m in re.finditer(r"G[KR].[DNE]G", s)} |
                  {m.start() for m in re.finditer(r"GKNDG|GKGDG", s)})
    print(f"  {key:12s} len={len(s):4d}  motifs: "
          + ", ".join(f"{s[h:h+5]}@{h+1}" for h in hits))

# the mouse TBCB CAP-Gly we established spans 134-231; use it as the query to find
# the homologous region of TBCE by sliding alignment
mtbcb = idx["mouse_Tbcb"]["sequence"]; q = mtbcb[133:231]
mtbce = idx["mouse_Tbce"]["sequence"]
print(f"\nquery = mouse TBCB CAP-Gly 134-231 ({len(q)} aa)")
best = None
for start in range(0, len(mtbce)-len(q)+1, 1):
    seg = mtbce[start:start+len(q)+20]
    A, B = nw_affine(q, seg)
    ident = sum(1 for x,y in zip(A,B) if x == y and x != "-")
    if best is None or ident > best[0]: best = (ident, start, seg)
ident, start, seg = best
print(f"best local match in mouse TBCE at residue {start+1}: {ident}/{len(q)} identities "
      f"({100*ident/len(q):.1f}%)")
print(f"  TBCE {start+1}-{start+len(q)+20}: {seg}")

# --- is it resolved in 9M1M chain E?
f = pdbx.CIFFile.read("data/pdb/9M1M.cif")
arr = pdbx.get_structure(f, model=1, use_author_fields=True)
E = arr[(arr.chain_id == "E") & struc.filter_amino_acids(arr)]
ca = E[E.atom_name == "CA"]; o = np.argsort(ca.res_id)
eseq = "".join(AA3.get(n,"X") for n in ca.res_name[o]); erid = ca.res_id[o].astype(int)
print(f"\n9M1M chain E (TBCE): {len(eseq)} residues resolved, auth {erid[0]}-{erid[-1]}")
hits = sorted({m.start() for m in re.finditer(r"G[KR].[DNE]G", eseq)})
for h in hits:
    print(f"  motif '{eseq[h:h+5]}' at auth {erid[h]}-{erid[h+4]}")
# which part of chain E is nearest the alpha-tubulin C-terminus (last resolved residue 437)?
A = arr[(arr.chain_id == "a") & struc.filter_amino_acids(arr)]
cterm = A[A.res_id >= 430]
d = np.linalg.norm(E.coord[:,None,:]-cterm.coord[None,:,:], axis=-1)
near = np.where(d.min(axis=1) < 12.0)[0]
nres = sorted(set(E.res_id[near].tolist()))
print(f"\nTBCE residues within 12 A of alpha-tubulin 430-437 (the start of the tail): {nres}")
json.dump(dict(tbcb_capgly=[134,231], tbce_capgly_match_start=int(start+1),
               tbce_capgly_identity=int(ident),
               tbce_near_alpha_cterm=nres),
          open("results/structure/tbce_capgly.json","w"), indent=1)
