#!/usr/bin/env python3
"""Where do the 45 TUBA8-divergent positions sit, structurally?

For alpha-tubulin we compute, from experimental coordinates:
  relSASA_mono      - relative SASA of the residue in the isolated alpha-tubulin chain
  dSASA_beta        - SASA buried by the alpha-beta intradimer interface
  dSASA_cofactor    - SASA buried by TBCD/TBCE/TBCC/Arl2 in 9M1M
  d_nucleotide      - min heavy-atom distance to the alpha-site nucleotide / Mg
  d_<cofactor>      - min heavy-atom distance to each cofactor chain
and then test, by permutation over all resolved positions, whether the divergent set is
enriched at (a) the cofactor interface, (b) the alpha-beta interface, (c) the buried core.
"""
import json, csv, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
MAXSASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
           'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}

def res_sasa(arr):
    """per-residue SASA summed from atom SASA (Shrake-Rupley), keyed by res_id."""
    sasa = struc.sasa(arr, vdw_radii="Single", point_number=200)
    out = {}
    for rid, s in zip(arr.res_id, sasa):
        if np.isnan(s): continue
        out[int(rid)] = out.get(int(rid), 0.0) + float(s)
    return out

def main():
    f = pdbx.CIFFile.read("data/pdb/9M1M.cif")
    arr = pdbx.get_structure(f, model=1, use_author_fields=True)
    arr = arr[~np.isin(arr.res_name, ["HOH"])]
    aa = struc.filter_amino_acids(arr)

    A   = arr[(arr.chain_id == "a") & aa]                       # alpha-tubulin
    B   = arr[(arr.chain_id == "b") & aa]                       # beta-tubulin
    COF = {c: arr[(arr.chain_id == c) & aa] for c in ["D","E","C","G"]}
    COFNAME = {"D":"TBCD","E":"TBCE","C":"TBCC","G":"ARL2"}

    het = arr[~aa & ~np.isin(arr.res_name, ["HOH"])]
    print("hetero present:", sorted(set(het.res_name.tolist())))
    # alpha-site nucleotide = the non-protein group closest to chain a
    nuc = None
    for rn in sorted(set(het.res_name.tolist())):
        sel = het[het.res_name == rn]
        for rid in sorted(set(sel.res_id.tolist())):
            grp = sel[sel.res_id == rid]
            d = np.min(np.linalg.norm(A.coord[:,None,:]-grp.coord[None,:,:], axis=-1))
            if d < 4.0:
                print(f"  group {rn}{rid} within {d:.2f} A of alpha-tubulin")
                nuc = grp if nuc is None else nuc + grp

    sasa_mono = res_sasa(A)
    sasa_ab   = res_sasa(A + B)                       # SASA of a in presence of b (take a's residues)
    sasa_full = res_sasa(A + B + sum(list(COF.values())[1:], list(COF.values())[0]))
    # note: res_sasa keys by res_id; chains a/b/cofactors share res_id ranges, so recompute carefully
    def res_sasa_of_chain(context, chain):
        s = struc.sasa(context, vdw_radii="Single", point_number=200)
        out = {}
        for ch, rid, v in zip(context.chain_id, context.res_id, s):
            if ch != chain or np.isnan(v): continue
            out[int(rid)] = out.get(int(rid), 0.0) + float(v)
        return out
    cof_all = None
    for c in COF.values(): cof_all = c if cof_all is None else cof_all + c
    sasa_mono = res_sasa_of_chain(A, "a")
    sasa_ab   = res_sasa_of_chain(A + B, "a")
    sasa_full = res_sasa_of_chain(A + B + cof_all, "a")

    aa1 = {}
    for rid, rn in zip(A.res_id, A.res_name): aa1[int(rid)] = AA3.get(rn, "X")

    rows = []
    for rid in sorted(sasa_mono):
        aa_i = aa1.get(rid, "X")
        mx = MAXSASA.get(aa_i, 200.0)
        sub = A[A.res_id == rid]
        rec = {"canon_pos": rid, "aa_TUBA1B": aa_i,
               "relSASA_mono": round(sasa_mono[rid]/mx, 3),
               "dSASA_beta": round(sasa_mono[rid]-sasa_ab.get(rid, 0.0), 2),
               "dSASA_cofactor": round(sasa_ab.get(rid,0.0)-sasa_full.get(rid,0.0), 2)}
        if nuc is not None:
            rec["d_nucleotide"] = round(float(np.min(np.linalg.norm(
                sub.coord[:,None,:]-nuc.coord[None,:,:], axis=-1))), 2)
        for c, cs in COF.items():
            rec[f"d_{COFNAME[c]}"] = round(float(np.min(np.linalg.norm(
                sub.coord[:,None,:]-cs.coord[None,:,:], axis=-1))), 2)
        rec["d_beta"] = round(float(np.min(np.linalg.norm(
            sub.coord[:,None,:]-B.coord[None,:,:], axis=-1))), 2)
        rows.append(rec)

    os.makedirs("results/structure", exist_ok=True)
    with open("results/structure/alpha_residue_structural_context_9M1M.csv","w",newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    json.dump(rows, open("results/structure/alpha_residue_structural_context_9M1M.json","w"), indent=1)
    print(f"\nwrote structural context for {len(rows)} resolved alpha-tubulin residues")
    return rows

if __name__ == "__main__":
    main()
