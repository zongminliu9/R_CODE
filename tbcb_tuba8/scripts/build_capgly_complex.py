#!/usr/bin/env python3
"""Transplant the experimental EEY binding mode (2HL3 chain C) into the mouse TBCB CAP-Gly
groove (1WHG) and, as an independent check, into the high-resolution T. brucei TBCB CAP-Gly
(4B6M, 1.59 A).

Superposition is by iterative Kabsch on sequence-aligned CA pairs with outlier rejection, so a
few divergent loops cannot drag the fit. The GKNDG motif -- the groove wall that actually reads
the -EEY/F-COO- terminus, and which is strictly conserved between p150Glued and mouse TBCB --
is used to verify the result rather than to drive it.

1WHG is an NESG construct: GSSGSSG at the N terminus and SGPSSG at the C terminus are tags and
are removed. The retained core is mapped onto full-length mouse TBCB (UniProt Q9D1E6) so that
residue numbers reported later refer to the real protein.
"""
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdbio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_isoforms import nw_affine

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
TAGS = ("GSSGSSG", "SGPSSG", "GDSGPSSG")

def load(cif):
    f = pdbx.CIFFile.read(cif)
    return pdbx.get_structure(f, model=1, use_author_fields=True)

def seq_of(arr, chain):
    sel = arr[(arr.chain_id == chain) & struc.filter_amino_acids(arr)]
    ca = sel[sel.atom_name == "CA"]; o = np.argsort(ca.res_id)
    return "".join(AA3.get(n,"X") for n in ca.res_name[o]), ca.res_id[o].astype(int), ca.coord[o].astype(float)

def kabsch(P, Q):
    """rotation+translation mapping P onto Q (both (n,3))"""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P-Pc).T @ (Q-Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, Qc - R @ Pc

def superpose_iter(P, Q, n_iter=12, cut0=6.0):
    """iterative Kabsch with outlier rejection; returns R,t,rmsd,mask"""
    mask = np.ones(len(P), bool); cut = cut0
    for it in range(n_iter):
        R, t = kabsch(P[mask], Q[mask])
        d = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
        rmsd = np.sqrt((d[mask]**2).mean())
        new = d < max(cut, 1.5)
        if new.sum() < 12: break
        if np.array_equal(new, mask) and it > 2: break
        mask = new; cut = max(cut*0.82, 2.0)
    R, t = kabsch(P[mask], Q[mask])
    d = np.linalg.norm((R @ P.T).T + t - Q, axis=1)
    return R, t, float(np.sqrt((d[mask]**2).mean())), mask

def aligned_pairs(sA, ridA, sB, ridB):
    A, B = nw_affine(sA, sB); i = j = 0; out = []
    for x, y in zip(A, B):
        if x != "-": i += 1
        if y != "-": j += 1
        if x != "-" and y != "-": out.append((i-1, j-1, x, y))
    return out

os.makedirs("data/models", exist_ok=True)
ref = load("data/pdb/2HL3.cif")
sR, ridR, caR = seq_of(ref, "A")
pep = ref[(ref.chain_id == "C") & struc.filter_amino_acids(ref)]
print(f"2HL3 reference: CAP-Gly {len(sR)} res; peptide chain C = "
      f"{''.join(AA3.get(n,'X') for n in sorted(set(zip(pep.res_id,pep.res_name)))[0:1])}"
      f"{''.join(AA3.get(n,'X') for _,n in sorted(set(zip(pep.res_id,pep.res_name))))} "
      f"({pep.array_length()} atoms, res {pep.res_id.min()}-{pep.res_id.max()})")

# full-length mouse TBCB for real numbering
idx = json.load(open("data/meta/sequence_index.json"))
mtbcb = idx["mouse_Tbcb"]["sequence"]
print(f"mouse TBCB (Q9D1E6) length {len(mtbcb)}")

results = {}
for pid, chain, label in [("1WHG","A","mouse TBCB CAP-Gly (NMR)"),
                          ("4B6M","A","T. brucei TBCB CAP-Gly 1.59A")]:
    arr = load(f"data/pdb/{pid}.cif")
    s, rid, ca = seq_of(arr, chain)
    # strip expression tags
    keep = np.ones(len(s), bool)
    for tg in TAGS:
        p = s.find(tg)
        while p >= 0:
            if p < 10 or p+len(tg) > len(s)-10:
                keep[p:p+len(tg)] = False
            p = s.find(tg, p+1)
    core_idx = np.where(keep)[0]
    s_core = "".join(s[i] for i in core_idx)
    rid_core = rid[core_idx]
    print(f"\n--- {pid} {label}")
    print(f"    {len(s)} res -> core {len(s_core)} res after tag removal (auth {rid_core[0]}-{rid_core[-1]})")
    if pid == "1WHG":
        off = mtbcb.find(s_core[:25])
        print(f"    core maps to full-length mouse TBCB at residue {off+1} "
              f"(so auth {rid_core[0]} == Q9D1E6 residue {off+1}); "
              f"CAP-Gly spans {off+1}-{off+len(s_core)}")
        results["mouse_tbcb_capgly_range"] = [off+1, off+len(s_core)]

    pairs = aligned_pairs(s_core, rid_core, sR, ridR)
    ident = sum(1 for _,_,x,y in pairs if x == y)
    P = np.array([ca[core_idx[i]] for i,_,_,_ in pairs])
    Q = np.array([caR[j]          for _,j,_,_ in pairs])
    R, t, rmsd, mask = superpose_iter(P, Q)
    print(f"    aligned {len(pairs)} CA, {ident} identical ({100*ident/len(pairs):.1f}%)")
    print(f"    iterative Kabsch: core RMSD {rmsd:.2f} A over {mask.sum()}/{len(pairs)} positions")

    # transform the EEY peptide from the 2HL3 frame into THIS domain's frame (inverse map)
    Rinv = R.T; tinv = -R.T @ t
    pep_new = pep.copy()
    pep_new.coord = (Rinv @ pep.coord.T).T + tinv
    pep_new.chain_id = np.full(pep_new.array_length(), "P")

    dom = arr[(arr.chain_id == chain) & struc.filter_amino_acids(arr)]
    dom = dom[np.isin(dom.res_id, rid_core)]
    # clash / contact check
    d = np.linalg.norm(dom.coord[:,None,:]-pep_new.coord[None,:,:], axis=-1)
    nclash = int((d < 2.0).sum()); ncont = int((d < 4.5).sum())
    contact_res = sorted(set(dom.res_id[np.where(d.min(axis=1) < 4.5)[0]].tolist()))
    print(f"    transplanted peptide: {nclash} atom pairs <2.0 A (clashes), "
          f"{ncont} pairs <4.5 A (contacts)")
    print(f"    contacting domain residues (auth): {contact_res}")
    gk = s_core.find("GKNDG"); gk = s_core.find("GKGDG") if gk < 0 else gk
    if gk >= 0:
        motif_rid = rid_core[gk:gk+5].tolist()
        inmotif = [r for r in contact_res if r in motif_rid]
        print(f"    GKNDG-type motif at auth {motif_rid} -> {len(inmotif)}/5 motif residues "
              f"contact the peptide: {inmotif}")
    merged = dom + pep_new
    pf = pdbio.PDBFile(); pf.set_structure(merged)
    out = f"data/models/{pid}_capgly_EEY.pdb"; pf.write(out)
    print(f"    wrote {out}")
    results[pid] = dict(label=label, rmsd=rmsd, n_aligned=len(pairs), n_core=int(mask.sum()),
                        pct_id=round(100*ident/len(pairs),1), clashes=nclash, contacts=ncont,
                        contact_residues=contact_res, model=out)
json.dump(results, open("results/structure/capgly_transplant.json","w"), indent=1)
