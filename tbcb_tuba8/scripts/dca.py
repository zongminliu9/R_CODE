#!/usr/bin/env python3
"""Mean-field Direct Coupling Analysis of alpha-tubulin from the 50000-sequence AF3 MSA.

Why this and not the elastic network. An ANM is harmonic and purely geometric: it reports how a
spring network built from one structure transmits force. It is blind to entropy-driven coupling,
to side-chain-mediated communication, and to any coupling that is functional rather than
mechanical. A negative ANM result therefore cannot exclude long-range coupling -- it only says
the harmonic mechanical channel is unremarkable.

DCA probes something categorically different: which residue pairs have co-varied over evolution
MORE than can be explained by transitive chains of other correlations. Inverting the covariance
matrix is what removes those indirect chains, which is the whole point -- raw correlation is
dominated by transitivity and would manufacture coupling everywhere.

Pipeline
  1. a3m -> match-state alignment (lowercase insertions dropped), gappy sequences removed
  2. sequence reweighting at 80% identity, via one-hot Gram matrix, giving Meff
  3. reweighted single-site and pair frequencies with pseudocount
  4. C = f_ij - f_i f_j over the 20 non-gap states, size (20L)^2
  5. e_ij = -(C^-1) block; score = Frobenius norm with average-product correction (APC)
"""
import os, sys, json, time, warnings
warnings.filterwarnings("ignore")
import numpy as np

AAs = "ACDEFGHIKLMNPQRSTVWY"          # 20 states; gap handled as index 20
A2I = {a: i for i, a in enumerate(AAs)}
Q = 21
MSA = "/home/user/uploads_fold/a98a190e-fold_tuba1a/msas/fold_tuba1a_paired_msa_chains_a.a3m"
MAXSEQ = int(os.environ.get("DCA_MAXSEQ", "12000"))
SEED = 0

def read_a3m(path, maxseq):
    """Return (M,L) int8 array of match-state columns. a3m lowercase = insertion, dropped."""
    seqs = []; name = None; buf = []
    qlen = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if buf:
                    s = "".join(buf)
                    m = "".join(c for c in s if not c.islower())
                    if qlen is None: qlen = len(m)
                    if len(m) == qlen: seqs.append(m)
                buf = []; name = line
            else:
                buf.append(line)
        if buf:
            s = "".join(buf); m = "".join(c for c in s if not c.islower())
            if qlen is not None and len(m) == qlen: seqs.append(m)
    L = qlen
    print(f"  parsed {len(seqs)} match-state sequences of length {L}")
    arr = np.full((len(seqs), L), Q-1, dtype=np.int8)
    for i, s in enumerate(seqs):
        for j, c in enumerate(s):
            arr[i, j] = A2I.get(c, Q-1)
    # drop sequences with >50% gaps
    gapfrac = (arr == Q-1).mean(axis=1)
    keep = gapfrac < 0.5
    arr = arr[keep]
    print(f"  after gap filter (<50% gaps): {arr.shape[0]} sequences")
    if arr.shape[0] > maxseq:
        rng = np.random.default_rng(SEED)
        # always keep the query (row 0)
        sel = np.concatenate([[0], rng.choice(np.arange(1, arr.shape[0]),
                                              size=maxseq-1, replace=False)])
        arr = arr[sel]
        print(f"  subsampled to {arr.shape[0]} sequences (query retained)")
    return arr, L

def onehot(arr, L):
    M = arr.shape[0]
    X = np.zeros((M, L*Q), dtype=np.float32)
    rows = np.repeat(np.arange(M), L)
    cols = (np.tile(np.arange(L), M)*Q + arr.reshape(-1)).astype(np.int64)
    X[rows, cols] = 1.0
    return X

def weights(X, L, theta=0.8):
    """w_i = 1/|{j: identity(i,j) >= theta}| from the one-hot Gram matrix."""
    t0 = time.time()
    S = X @ X.T                       # identity counts (M,M)
    thr = theta*L
    nsim = (S >= thr).sum(axis=1).astype(np.float64)
    w = 1.0/np.maximum(nsim, 1.0)
    print(f"  reweighting done in {time.time()-t0:.1f}s; Meff = {w.sum():.1f} "
          f"(M = {X.shape[0]}, Meff/L = {w.sum()/L:.2f})")
    return w

def freqs(arr, w, L, pseudo=0.5):
    Meff = w.sum()
    lam = pseudo*Meff
    X = onehot(arr, L)
    Wm = (w[:, None].astype(np.float32))
    fi = (X*Wm).sum(axis=0)/ (Meff+lam)
    fi = fi.reshape(L, Q) + lam/((Meff+lam)*Q)
    t0 = time.time()
    Xw = X*np.sqrt(Wm)
    Fij = (Xw.T @ Xw)/(Meff+lam)
    print(f"  pair frequencies in {time.time()-t0:.1f}s  ({Fij.shape})")
    Fij = Fij.reshape(L, Q, L, Q) + lam/((Meff+lam)*Q*Q)
    # restore exact single-site marginals on the diagonal blocks
    for i in range(L):
        Fij[i,:,i,:] = np.diag(fi[i])
    return fi, Fij

def mfdca(fi, Fij, L):
    q1 = Q-1
    C = np.empty((L*q1, L*q1), dtype=np.float64)
    for i in range(L):
        for j in range(L):
            C[i*q1:(i+1)*q1, j*q1:(j+1)*q1] = (Fij[i,:q1,j,:q1]
                                               - np.outer(fi[i,:q1], fi[j,:q1]))
    t0 = time.time()
    P = np.linalg.inv(C)
    print(f"  inverted {C.shape[0]}x{C.shape[0]} covariance in {time.time()-t0:.1f}s")
    FN = np.zeros((L, L))
    for i in range(L):
        bi = slice(i*q1, (i+1)*q1)
        for j in range(i+1, L):
            e = -P[bi, j*q1:(j+1)*q1]
            # zero-sum (Ising) gauge
            e = e - e.mean(axis=0, keepdims=True) - e.mean(axis=1, keepdims=True) + e.mean()
            v = np.sqrt((e**2).sum())
            FN[i, j] = FN[j, i] = v
    # average-product correction removes per-site entropic bias
    s = FN.sum(axis=0); tot = FN.sum()
    APC = np.outer(s, s)/max(tot, 1e-30)
    DI = FN - APC
    np.fill_diagonal(DI, 0.0)
    return FN, DI

if __name__ == "__main__":
    os.makedirs("results/dca", exist_ok=True)
    print("reading MSA ...")
    arr, L = read_a3m(MSA, MAXSEQ)
    X = onehot(arr, L)
    w = weights(X, L)
    del X
    fi, Fij = freqs(arr, w, L)
    FN, DI = mfdca(fi, Fij, L)
    np.save("results/dca/dca_FN.npy", FN.astype(np.float32))
    np.save("results/dca/dca_APC.npy", DI.astype(np.float32))
    json.dump(dict(L=int(L), M=int(arr.shape[0]), Meff=float(w.sum()),
                   Meff_over_L=float(w.sum()/L), msa=MSA, maxseq=MAXSEQ),
              open("results/dca/dca_meta.json","w"), indent=1)
    # sanity: top APC-corrected pairs should be mostly spatial contacts
    iu = np.triu_indices(L, k=5)
    vals = DI[iu]; order = np.argsort(vals)[::-1][:25]
    print("\n  top 25 APC-corrected couplings (|i-j|>=5):")
    for k in order[:25]:
        print(f"    {iu[0][k]+1:4d} - {iu[1][k]+1:4d}   {vals[k]:.4f}")
