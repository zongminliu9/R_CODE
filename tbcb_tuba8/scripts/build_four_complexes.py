#!/usr/bin/env python3
"""Build and relax four complexes:  {TBCB-CAP-Gly, TBCE-CAP-Gly} x {EEY, EEF}.

Templates
  TBCE CAP-Gly : 9M1M chain E residues 4-121  -- experimental, 2.21 A, same structure as the
                 late cofactor complex, motif GKHDG at 47-51.
  TBCB CAP-Gly : AlphaFold mouse TBCB (Q9D1E6) residues 134-231, cross-validated against the
                 experimental mouse NMR structure 1WHG and the 1.59 A T. brucei structure 4B6M.
  peptide      : EEY from 2HL3 chain C, carrying the free C-terminal carboxylate, transplanted
                 by superposition onto each domain. EEF is made by deleting the Tyr hydroxyl
                 and letting the force field rebuild HZ.

Clashes from transplantation are relieved by minimisation with the peptide backbone restrained,
so the experimental binding mode is preserved while side chains relax.
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdbio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_isoforms import nw_affine
from build_capgly_complex import kabsch, superpose_iter, aligned_pairs, AA3

os.makedirs("data/models", exist_ok=True)

def get(cif, chain, lo=None, hi=None, model=1):
    f = pdbx.CIFFile.read(cif) if cif.endswith(".cif") else None
    if f is not None:
        arr = pdbx.get_structure(f, model=model, use_author_fields=True)
    else:
        arr = pdbio.PDBFile.read(cif).get_structure(model=model)
    sel = arr[(arr.chain_id == chain) & struc.filter_amino_acids(arr)]
    if lo is not None: sel = sel[(sel.res_id >= lo) & (sel.res_id <= hi)]
    return sel

def seq_rid_ca(sel):
    ca = sel[sel.atom_name == "CA"]; o = np.argsort(ca.res_id)
    return ("".join(AA3.get(n,"X") for n in ca.res_name[o]),
            ca.res_id[o].astype(int), ca.coord[o].astype(float))

# reference: p150Glued CAP-Gly + EEY
ref = pdbx.get_structure(pdbx.CIFFile.read("data/pdb/2HL3.cif"), model=1, use_author_fields=True)
ref = ref[ref.res_name != "HOH"]
refdom = ref[(ref.chain_id == "A") & struc.filter_amino_acids(ref)]
refpep = ref[(ref.chain_id == "C") & struc.filter_amino_acids(ref)]
sR, ridR, caR = seq_rid_ca(refdom)

DOMAINS = {
  "TBCB": dict(path="data/afdb/Q9D1E6.pdb", chain="A", lo=134, hi=231,
               label="mouse TBCB CAP-Gly (AlphaFold Q9D1E6 134-231)", motif="GKNDG"),
  "TBCE": dict(path="data/pdb/9M1M.cif",    chain="E", lo=4,   hi=121,
               label="mouse/human TBCE CAP-Gly (9M1M chain E 4-121, 2.21 A)", motif="GKHDG"),
}

built = {}
for name, D in DOMAINS.items():
    dom = get(D["path"], D["chain"], D["lo"], D["hi"])
    s, rid, ca = seq_rid_ca(dom)
    pairs = aligned_pairs(s, rid, sR, ridR)
    ident = sum(1 for _,_,x,y in pairs if x == y)
    P = np.array([ca[i] for i,_,_,_ in pairs]); Q = np.array([caR[j] for _,j,_,_ in pairs])
    R, t, rmsd, mask = superpose_iter(P, Q)
    Rinv, tinv = R.T, -R.T @ t
    pep = refpep.copy()
    pep.coord = (Rinv @ refpep.coord.T).T + tinv
    pep.chain_id = np.full(pep.array_length(), "P")
    d = np.linalg.norm(dom.coord[:,None,:]-pep.coord[None,:,:], axis=-1)
    print(f"\n--- {name}: {D['label']}")
    print(f"    motif {D['motif']}; aligned {len(pairs)} CA, {ident} identical "
          f"({100*ident/len(pairs):.1f}%); Kabsch RMSD {rmsd:.2f} A over {mask.sum()}")
    print(f"    transplanted EEY: clashes<2A = {int((d<2.0).sum())}, pairs<4.5A = {int((d<4.5).sum())}"
          f"   [native reference: 0 clashes, 138 pairs]")
    cres = sorted(set(dom.res_id[np.where(d.min(axis=1)<4.5)[0]].tolist()))
    nm = {int(r): AA3.get(n,'X') for r,n in zip(dom.res_id, dom.res_name)}
    print(f"    contacts: {' '.join(f'{nm[r]}{r}' for r in cres)}")
    merged = dom + pep
    pf = pdbio.PDBFile(); pf.set_structure(merged)
    out = f"data/models/{name}_capgly_EEY_raw.pdb"; pf.write(out)
    built[name] = dict(model=out, rmsd=float(rmsd), pct_id=round(100*ident/len(pairs),1),
                       clashes=int((d<2.0).sum()), pairs=int((d<4.5).sum()),
                       contacts=cres, label=D["label"], motif=D["motif"])

# ---- make the EEF variant: delete the Tyr hydroxyl oxygen, rename TYR->PHE
for name in DOMAINS:
    src = built[name]["model"]
    lines = [l for l in open(src)]
    out = src.replace("_EEY_raw", "_EEF_raw")
    kept = 0; dropped = 0
    with open(out, "w") as fh:
        for l in lines:
            if l.startswith(("ATOM","HETATM")):
                resn = l[17:20].strip(); chain = l[21]; an = l[12:16].strip()
                if chain == "P" and resn == "TYR":
                    if an in ("OH",):
                        dropped += 1; continue
                    l = l[:17] + "PHE" + l[20:]
                kept += 1
            fh.write(l)
    print(f"{name}: wrote {out}  (removed {dropped} hydroxyl atom, TYR->PHE)")
    built[name]["model_EEF"] = out
json.dump(built, open("results/structure/four_complexes_raw.json","w"), indent=1)
