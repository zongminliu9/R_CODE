#!/usr/bin/env python3
"""Ground-truth check: what does the EEY peptide actually touch in the 2HL3 crystal structure?
If the transplanted models do not reproduce that contact pattern, the transplant is misplaced
and must be fixed before any energetics are computed."""
import warnings, numpy as np, json
warnings.filterwarnings("ignore")
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx, biotite.structure.io.pdb as pdbio
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}

f = pdbx.CIFFile.read("data/pdb/2HL3.cif")
arr = pdbx.get_structure(f, model=1, use_author_fields=True)
arr = arr[arr.res_name != "HOH"]
aa = struc.filter_amino_acids(arr)
dom = arr[(arr.chain_id == "A") & aa]
pep = arr[(arr.chain_id == "C") & aa]
print(f"2HL3 native: domain chain A {len(set(dom.res_id))} res; "
      f"peptide chain C {[f'{AA3.get(n)}{r}' for r,n in sorted(set(zip(pep.res_id,pep.res_name)))]}")

def contacts(dom, pep, cut=4.5):
    d = np.linalg.norm(dom.coord[:,None,:]-pep.coord[None,:,:], axis=-1)
    hit = d.min(axis=1) < cut
    res = {}
    for i in np.where(hit)[0]:
        r = int(dom.res_id[i]); res.setdefault(r, []).append(float(d[i].min()))
    return {r: round(min(v),2) for r,v in sorted(res.items())}, int((d<2.0).sum()), int((d<cut).sum())

nat, ncl, nct = contacts(dom, pep)
names = {int(r): AA3.get(n,'X') for r,n in zip(dom.res_id, dom.res_name)}
print(f"\nNATIVE contact residues (<4.5 A), p150Glued numbering:")
print("   " + "  ".join(f"{names[r]}{r}({d})" for r,d in nat.items()))
print(f"   clashes<2.0A = {ncl}  (crystal structure, so this is the physical baseline)")
print(f"   total pairs<4.5A = {nct}")

# per-peptide-residue breakdown, and specifically what the terminal Tyr touches
print("\nwhat each peptide residue touches:")
for r in sorted(set(pep.res_id)):
    sub = pep[pep.res_id == r]
    d = np.linalg.norm(dom.coord[:,None,:]-sub.coord[None,:,:], axis=-1)
    part = sorted({int(dom.res_id[i]) for i in np.where(d.min(axis=1) < 4.5)[0]})
    nm = AA3.get(sub.res_name[0],'X')
    print(f"   {nm}{r}: {['%s%d'%(names[p],p) for p in part]}")
# the terminal carboxylate and the ring hydroxyl specifically
term = pep[pep.res_id == pep.res_id.max()]
for an in ("OXT","O","OH","CZ"):
    s = term[term.atom_name == an]
    if s.array_length() == 0: continue
    d = np.linalg.norm(dom.coord - s.coord[0], axis=1)
    near = np.argsort(d)[:6]
    print(f"   Tyr {an:3s} nearest domain atoms: "
          + ", ".join(f"{names[int(dom.res_id[i])]}{dom.res_id[i]}:{dom.atom_name[i]}={d[i]:.2f}" for i in near))

print("\n=== transplanted models vs native pattern ===")
aln = {}
for pid in ("1WHG","4B6M"):
    pf = pdbio.PDBFile.read(f"data/models/{pid}_capgly_EEY.pdb")
    m = pf.get_structure(model=1)
    md = m[m.chain_id != "P"]; mp = m[m.chain_id == "P"]
    c, cl, ct = contacts(md, mp)
    mn = {int(r): AA3.get(n,'X') for r,n in zip(md.res_id, md.res_name)}
    print(f"\n  {pid}: {len(c)} contact residues, clashes<2A = {cl}, pairs<4.5A = {ct}")
    print("     " + "  ".join(f"{mn.get(r,'X')}{r}" for r in c))
json.dump(dict(native_contacts=nat, native_clashes=ncl, native_pairs=nct),
          open("results/structure/native_2HL3_contacts.json","w"), indent=1)
