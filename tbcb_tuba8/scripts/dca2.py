#!/usr/bin/env python3
"""DCA, take two: all 44191 sequences via chunked reweighting, then VALIDATED against structure.

Two fixes over the first pass.

(1) Meff was only 515 (Meff/L = 1.14, well below the ~5 needed for mean-field DCA) because
    subsampling to 12000 threw away sequences and alpha-tubulin is so conserved that 80%-identity
    reweighting collapses near-duplicates. Using all 44191 requires the M x M identity matrix,
    which is 7.8 GB at M=44191, so it is computed in row blocks and never materialised.

(2) Nothing should be interpreted before the method is shown to work on this family. The standard
    test is contact precision: the top-ranked long-range couplings of a working DCA land on
    spatial contacts in the known structure. If precision is near the random baseline the
    couplings are noise and no allosteric claim can rest on them.
"""
import os, sys, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx

AAs = "ACDEFGHIKLMNPQRSTVWY"; A2I = {a:i for i,a in enumerate(AAs)}; Q = 21
MSA = "/home/user/uploads_fold/a98a190e-fold_tuba1a/msas/fold_tuba1a_paired_msa_chains_a.a3m"

def read_a3m(path):
    seqs=[]; buf=[]; qlen=None
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if buf:
                    m="".join(c for c in "".join(buf) if not c.islower())
                    if qlen is None: qlen=len(m)
                    if len(m)==qlen: seqs.append(m)
                buf=[]
            else: buf.append(line)
        if buf:
            m="".join(c for c in "".join(buf) if not c.islower())
            if qlen is not None and len(m)==qlen: seqs.append(m)
    L=qlen
    arr=np.full((len(seqs),L), Q-1, dtype=np.int8)
    for i,s in enumerate(seqs):
        for j,c in enumerate(s): arr[i,j]=A2I.get(c,Q-1)
    keep=(arr==Q-1).mean(axis=1) < 0.5
    arr=arr[keep]
    print(f"  {arr.shape[0]} sequences x {L} match-state columns")
    return arr, L

def onehot(arr, L, dtype=np.float32):
    M=arr.shape[0]; X=np.zeros((M,L*Q), dtype=dtype)
    X[np.repeat(np.arange(M),L), (np.tile(np.arange(L),M)*Q+arr.reshape(-1)).astype(np.int64)]=1.0
    return X

def weights_chunked(X, L, theta=0.8, block=1500):
    """1/|{j : seq-identity(i,j) >= theta}| without ever forming the full M x M matrix."""
    M=X.shape[0]; thr=theta*L; nsim=np.zeros(M)
    t0=time.time()
    for s in range(0, M, block):
        e=min(s+block, M)
        S = X[s:e] @ X.T
        nsim[s:e] = (S >= thr).sum(axis=1)
        del S
    w = 1.0/np.maximum(nsim, 1.0)
    print(f"  reweighting {time.time()-t0:.1f}s   Meff = {w.sum():.1f}   Meff/L = {w.sum()/L:.2f}")
    return w

def dca(arr, w, L, pseudo=0.5):
    Meff=w.sum(); lam=pseudo*Meff
    X=onehot(arr,L)
    Wm=w[:,None].astype(np.float32)
    fi=(X*Wm).sum(axis=0)/(Meff+lam); fi=fi.reshape(L,Q)+lam/((Meff+lam)*Q)
    Xw=X*np.sqrt(Wm); del X
    t0=time.time(); Fij=(Xw.T@Xw)/(Meff+lam); del Xw
    print(f"  pair frequencies {time.time()-t0:.1f}s")
    Fij=Fij.reshape(L,Q,L,Q)+lam/((Meff+lam)*Q*Q)
    for i in range(L): Fij[i,:,i,:]=np.diag(fi[i])
    q1=Q-1; C=np.empty((L*q1,L*q1))
    for i in range(L):
        for j in range(L):
            C[i*q1:(i+1)*q1, j*q1:(j+1)*q1]=Fij[i,:q1,j,:q1]-np.outer(fi[i,:q1],fi[j,:q1])
    del Fij
    t0=time.time(); P=np.linalg.inv(C); del C
    print(f"  matrix inversion {time.time()-t0:.1f}s")
    FN=np.zeros((L,L))
    for i in range(L):
        bi=slice(i*q1,(i+1)*q1)
        for j in range(i+1,L):
            e=-P[bi, j*q1:(j+1)*q1]
            e=e-e.mean(0,keepdims=True)-e.mean(1,keepdims=True)+e.mean()
            FN[i,j]=FN[j,i]=np.sqrt((e**2).sum())
    s=FN.sum(0); DI=FN-np.outer(s,s)/FN.sum(); np.fill_diagonal(DI,0.0)
    return FN, DI

def structure_contacts(cutoff=8.0):
    """Cb-Cb (Ca for Gly) contact map of alpha-tubulin from 9M1M chain a, canonical numbering."""
    arr=pdbx.get_structure(pdbx.CIFFile.read("data/pdb/9M1M.cif"), model=1, use_author_fields=True)
    A=arr[(arr.chain_id=="a") & struc.filter_amino_acids(arr)]
    pos={}
    for rid in np.unique(A.res_id):
        sub=A[A.res_id==rid]
        cb=sub[sub.atom_name=="CB"]
        ca=sub[sub.atom_name=="CA"]
        p = cb.coord[0] if cb.array_length() else (ca.coord[0] if ca.array_length() else None)
        if p is not None: pos[int(rid)]=p
    ks=sorted(pos); X=np.array([pos[k] for k in ks])
    D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=-1)
    return ks, D < cutoff, D

if __name__=="__main__":
    os.makedirs("results/dca", exist_ok=True)
    arr,L=read_a3m(MSA)
    X=onehot(arr,L); w=weights_chunked(X,L); del X
    FN,DI=dca(arr,w,L)
    np.save("results/dca/dca_FN_full.npy", FN.astype(np.float32))
    np.save("results/dca/dca_APC_full.npy", DI.astype(np.float32))
    meta=dict(L=int(L), M=int(arr.shape[0]), Meff=float(w.sum()), Meff_over_L=float(w.sum()/L))
    json.dump(meta, open("results/dca/dca_meta_full.json","w"), indent=1)

    print("\n=== VALIDATION: contact precision against 9M1M chain a (Cb-Cb < 8 A) ===")
    ks, cmap, D = structure_contacts()
    kset={k:i for i,k in enumerate(ks)}
    rows=[]
    for sep,label in ((6,"|i-j|>=6"),(12,"|i-j|>=12"),(24,"|i-j|>=24")):
        pairs=[]
        for i in range(L):
            for j in range(i+sep, L):
                p,q=i+1,j+1
                if p in kset and q in kset: pairs.append((DI[i,j], p, q))
        pairs.sort(reverse=True)
        base = np.mean([cmap[kset[p],kset[q]] for _,p,q in pairs])
        print(f"  {label}: {len(pairs)} scorable pairs, random baseline {100*base:5.2f}%")
        for n in (int(L*0.25), int(L*0.5), L, 2*L):
            top=pairs[:n]
            prec=np.mean([cmap[kset[p],kset[q]] for _,p,q in top])
            print(f"     top {n:4d}  precision {100*prec:6.2f}%   "
                  f"enrichment {prec/base:5.2f}x")
            rows.append(dict(sep=sep, n=n, precision=float(prec), baseline=float(base),
                             enrichment=float(prec/base)))
    json.dump(rows, open("results/dca/dca_validation.json","w"), indent=1)
