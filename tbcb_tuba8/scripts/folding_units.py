#!/usr/bin/env python3
"""Where does TUBA8's divergence sit in the FOLDING process, not the final state?

The objection that AF3 models an end state is right, and the way to act on it is not to guess at
intermediates but to ask which parts of alpha-tubulin fold as units and which fold late, then
locate TUBA8's substitutions on that map.

Three independent proxies, deliberately chosen because none of them needs a trajectory:

  1. CONTACT ORDER per residue -- the mean sequence separation of a residue's spatial contacts.
     Low contact order means a residue's environment is made of sequence-local partners, which
     forms early and independently; high contact order means it needs distant parts of the chain
     to be in place first, so it forms late. This is the same quantity whose chain average
     predicts folding rates across small proteins.

  2. FOLDING UNITS by contact-map modularity -- spectral clustering of the residue contact graph
     gives the compact, mostly self-contacting blocks that behave as cooperative folding units
     (foldons). A substitution inside a unit perturbs that unit's own folding; one at a
     unit-unit boundary perturbs their assembly, which is a later and more fragile step.

  3. COEVOLUTIONARY DEGREE from the DCA already computed -- coupling hubs tend to be the
     residues whose identity the rest of the fold depends on, which is the evolutionary
     signature of a folding nucleus.

The question each answers: are TUBA8's substitutions in early-forming local structure (harmless
to the pathway, consistent with the buffering already found), or in the late, cooperative steps
where a chaperone would be needed?
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
import biotite.structure as struc, biotite.structure.io.pdbx as pdbx

arr = pdbx.get_structure(pdbx.CIFFile.read("data/pdb/9M1M.cif"), model=1, use_author_fields=True)
A = arr[(arr.chain_id=="a") & struc.filter_amino_acids(arr)]
res = sorted(np.unique(A.res_id).tolist())
N = len(res); r2i = {r:i for i,r in enumerate(res)}
# heavy-atom contact map at 5 A (any-atom), which is the standard definition for contact order
coords = [A[A.res_id==r].coord.astype(float) for r in res]
C = np.zeros((N,N), bool)
for i in range(N):
    for j in range(i+1, N):
        if np.min(np.linalg.norm(coords[i][:,None,:]-coords[j][None,:,:], axis=-1)) < 5.0:
            C[i,j]=C[j,i]=True
print(f"alpha-tubulin chain a: {N} residues, {int(C.sum()/2)} heavy-atom contacts <5 A")

# ---------- 1. contact order
sep = np.abs(np.subtract.outer(np.array(res), np.array(res)))
co = np.array([ sep[i][C[i]].mean() if C[i].any() else 0.0 for i in range(N) ])
rco = co/N
print(f"chain-average relative contact order = {rco.mean():.4f} "
      f"(globular proteins of this size: 0.08-0.15)")

# ---------- 2. folding units by spectral clustering of the contact graph
W = C.astype(float) + 0.3*np.eye(N)*0        # pure contact graph
# add weak chain connectivity so units are contiguous-ish but not forced
for i in range(N-1): W[i,i+1]=W[i+1,i]=max(W[i,i+1],1.0)
deg = W.sum(1); Lap = np.diag(deg) - W
Dm = np.diag(1/np.sqrt(np.maximum(deg,1e-9)))
Lsym = Dm @ Lap @ Dm
vals, vecs = np.linalg.eigh(Lsym)
K = 6
emb = vecs[:, 1:K]
emb = emb/np.maximum(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12)
from scipy.cluster.vq import kmeans2
cent, lab = kmeans2(emb, K, minit='++', seed=0)
print(f"\nspectral folding units (K={K}):")
for k in range(K):
    mem = [res[i] for i in range(N) if lab[i]==k]
    if not mem: continue
    intra = C[np.ix_([r2i[m] for m in mem],[r2i[m] for m in mem])].sum()/2
    tot   = C[[r2i[m] for m in mem]].sum() - intra
    segs=[]; s=mem[0]; p=mem[0]
    for m in mem[1:]:
        if m==p+1: p=m; continue
        segs.append((s,p)); s=p=m
    segs.append((s,p))
    segs=[f"{a}-{b}" for a,b in segs if b-a>=3]
    print(f"  unit {k}: {len(mem):3d} res, intra/inter contacts {intra:5.0f}/{tot:5.0f} "
          f"({intra/max(intra+tot,1):.2f} self)  segments {segs[:6]}")

# boundary residues: contacts to a different unit
bnd = np.array([ (lab[C[i]] != lab[i]).sum()/max(C[i].sum(),1) for i in range(N) ])

# ---------- 3. coevolutionary degree
DI = np.load("results/dca/dca_APC_full.npy").astype(np.float64)
deg_dca = np.array([ np.sort(DI[r-1])[-10:].sum() for r in res ])   # top-10 coupling strength

# ---------- locate TUBA8's divergence on all three
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
DIV = [d["canon_pos"] for d in div if d["canon_pos"] in r2i]
LOOP = [p for p in range(32,46) if p in r2i]
idx = np.array([r2i[p] for p in DIV]); lidx = np.array([r2i[p] for p in LOOP])
rng = np.random.default_rng(11)

print(f"\n{'metric':28s} {'divergent':>10s} {'rest':>8s} {'p(high)':>8s} {'p(low)':>8s}")
for nm, v in (("relative contact order", rco), ("unit-boundary fraction", bnd),
              ("coevolutionary degree", deg_dca)):
    obs = v[idx].mean(); rest = np.delete(v, idx).mean()
    null = np.array([v[rng.choice(N, size=len(idx), replace=False)].mean() for _ in range(50000)])
    print(f"{nm:28s} {obs:10.4f} {rest:8.4f} "
          f"{(np.sum(null>=obs)+1)/50001:8.4f} {(np.sum(null<=obs)+1)/50001:8.4f}")

print(f"\n{'metric':28s} {'32-45 loop':>11s} {'rest':>8s} {'p(high)':>8s} {'p(low)':>8s}")
for nm, v in (("relative contact order", rco), ("unit-boundary fraction", bnd),
              ("coevolutionary degree", deg_dca)):
    obs = v[lidx].mean()
    null = np.array([v[rng.choice(N, size=len(lidx), replace=False)].mean() for _ in range(50000)])
    print(f"{nm:28s} {obs:11.4f} {np.delete(v,lidx).mean():8.4f} "
          f"{(np.sum(null>=obs)+1)/50001:8.4f} {(np.sum(null<=obs)+1)/50001:8.4f}")

print("\nwhich folding unit carries the divergence?")
for k in range(K):
    mem = [res[i] for i in range(N) if lab[i]==k]
    nd = len(set(mem) & set(DIV))
    print(f"  unit {k}: {len(mem):3d} residues, {nd:2d} divergent "
          f"({100*nd/max(len(mem),1):4.1f}%)  loop32-45 members: "
          f"{len(set(mem)&set(LOOP))}")
json.dump(dict(residues=res, rel_contact_order=rco.tolist(), unit=lab.tolist(),
               boundary_frac=bnd.tolist(), dca_degree=deg_dca.tolist()),
          open("results/structure/folding_units.json","w"), indent=1)
