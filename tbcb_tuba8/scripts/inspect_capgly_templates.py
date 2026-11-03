#!/usr/bin/env python3
"""Report polymer entities (chain, length, one-letter seq) + ligands for candidate CAP-Gly complexes."""
import biotite.structure.io.pdbx as pdbx
import warnings; warnings.filterwarnings("ignore")

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'m','SEP':'s','TPO':'t','PTR':'y'}

for pid in ["2HL3","2HQH","2HKN","9YNE","9YNH","8PR0","3RDV","1WHG","4B6M","1V6E"]:
    try:
        f = pdbx.CIFFile.read(f"data/pdb/{pid}.cif")
        blk = f.block
        try:   title = blk["struct"]["title"].as_item()
        except Exception: title = ""
        arr = pdbx.get_structure(f, model=1, use_author_fields=True)
        print(f"\n--- {pid}  {str(title)[:76]}")
        import numpy as np
        for ch in np.unique(arr.chain_id):
            sub = arr[arr.chain_id == ch]
            # unique residues in order
            seen = {}; order = []
            for rid, rn, het in zip(sub.res_id, sub.res_name, sub.hetero):
                k = (rid, rn)
                if k not in seen: seen[k] = het; order.append(k)
            poly = [(r, n) for (r, n) in order if n in AA3]
            lig  = sorted({n for (r, n) in order if n not in AA3 and n not in ('HOH',)})
            if poly:
                seq = "".join(AA3[n] for _, n in poly)
                rng = f"{poly[0][0]}-{poly[-1][0]}"
                tag = "  <<< PEPTIDE" if len(poly) <= 15 else ""
                print(f"    chain {ch}: {len(poly):4d} res  [{rng:>10s}]  {seq[:70]}{'...' if len(seq)>70 else ''}{tag}")
                if len(poly) <= 15: print(f"                 full peptide seq = {seq}   C-term residue = {poly[-1][1]}")
            if lig: print(f"    chain {ch}: ligands {lig}")
    except Exception as e:
        print(f"\n--- {pid}  ERROR {type(e).__name__}: {e}")
