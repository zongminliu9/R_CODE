#!/usr/bin/env python3
"""Anisotropic network model + perturbation response scanning on alpha-tubulin.

Honest scope statement first. An ANM is built from Calpha geometry alone, so it is
sequence-agnostic: it CANNOT tell us that TUBA8 is intrinsically more mobile than TUBA1A,
because 90%-identical models have near-identical contact topology and any difference would
mostly reflect model noise. What it CAN do rigorously, from one high-resolution structure, is
answer a different question: do the positions where TUBA8 differs occupy mechanically special
sites in the common alpha-tubulin fold, and is the 32-45 loop mechanically coupled to the
functional surfaces it does not physically touch?

Method
  Hessian  H_ij = -gamma * (r_ij r_ij^T)/|r_ij|^2  for |r_ij| < cutoff, diagonal by row sum.
  Covariance C = pinv(H) over the 3N-6 non-trivial modes.
  MSF_i = Tr(C_ii).
  PRS: averaging unit random forces at i over directions gives
       S_ij = <|dR_j|^2> = (1/3) * ||C_ij||_F^2 .
  effectiveness(i) = mean_j S_ij   (how far i's push propagates)
  sensitivity(j)   = mean_i S_ij   (how much j reports others' pushes)
"""
import json, os, csv, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

CUTOFF = 13.0
GAMMA  = 1.0
rng = np.random.default_rng(11)

def load_ca(cif, chain):
    f = pdbx.CIFFile.read(cif)
    arr = pdbx.get_structure(f, model=1, use_author_fields=True)
    sel = arr[(arr.chain_id == chain) & (arr.atom_name == "CA") & struc.filter_amino_acids(arr)]
    order = np.argsort(sel.res_id)
    return sel.coord[order].astype(float), sel.res_id[order].astype(int)

def anm(coords, cutoff=CUTOFF, gamma=GAMMA):
    N = len(coords)
    H = np.zeros((3*N, 3*N))
    d = coords[:, None, :] - coords[None, :, :]
    dist = np.linalg.norm(d, axis=-1)
    for i in range(N):
        for j in range(i+1, N):
            r = dist[i, j]
            if r >= cutoff or r < 1e-8: continue
            v = d[i, j]/r
            k = gamma*np.outer(v, v)
            H[3*i:3*i+3, 3*j:3*j+3] -= k
            H[3*j:3*j+3, 3*i:3*i+3] -= k
            H[3*i:3*i+3, 3*i:3*i+3] += k
            H[3*j:3*j+3, 3*j:3*j+3] += k
    return H

def covariance(H, n_trivial=6):
    w, V = np.linalg.eigh(H)
    idx = np.argsort(w)
    w, V = w[idx], V[:, idx]
    keep = slice(n_trivial, len(w))
    wk, Vk = w[keep], V[:, keep]
    good = wk > 1e-8
    C = (Vk[:, good]/wk[good]) @ Vk[:, good].T
    return C, w, V

def prs(C, N):
    S = np.zeros((N, N))
    for i in range(N):
        Ci = C[:, 3*i:3*i+3]                       # (3N, 3)
        blocks = Ci.reshape(N, 3, 3)
        S[i] = (blocks**2).sum(axis=(1, 2))/3.0
    return S

if __name__ == "__main__":
    os.makedirs("results/anm", exist_ok=True)
    coords, resid = load_ca("data/pdb/9M1M.cif", "a")
    N = len(coords)
    print(f"alpha-tubulin chain a: N={N} CA atoms, auth {resid[0]}-{resid[-1]}, cutoff {CUTOFF} A")
    H = anm(coords)
    C, w, V = covariance(H)
    nz = w[6:][w[6:] > 1e-8]
    print(f"Hessian {H.shape}; first 6 eigenvalues (should be ~0): {np.array2string(w[:6], precision=1)}")
    print(f"first 6 non-trivial eigenvalues: {np.array2string(nz[:6], precision=4)}")
    msf = np.array([np.trace(C[3*i:3*i+3, 3*i:3*i+3]) for i in range(N)])
    S = prs(C, N)
    eff = S.mean(axis=1); sen = S.mean(axis=0)
    # normalise to z-scores for interpretability
    z = lambda x: (x-x.mean())/x.std()
    out = [dict(canon_pos=int(resid[i]), msf=float(msf[i]), msf_z=float(z(msf)[i]),
                effectiveness=float(eff[i]), eff_z=float(z(eff)[i]),
                sensitivity=float(sen[i]), sen_z=float(z(sen)[i])) for i in range(N)]
    with open("results/anm/alpha_anm_prs.csv","w",newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0].keys())); wr.writeheader(); wr.writerows(out)
    np.save("results/anm/prs_matrix.npy", S)
    np.save("results/anm/anm_eigenvalues.npy", w)
    json.dump(out, open("results/anm/alpha_anm_prs.json","w"), indent=1)
    print(f"wrote per-residue MSF / effectiveness / sensitivity for {N} residues")
    # quick look at the extremes
    pos = np.array([r["canon_pos"] for r in out])
    for nm, arr_ in (("effectiveness", eff), ("MSF", msf)):
        top = pos[np.argsort(arr_)[::-1][:12]]
        print(f"  highest {nm:13s}: {sorted(top.tolist())}")
