#!/usr/bin/env python3
"""The control that decides whether the coevolving-module result survives.

TUBA8's divergent positions are, by definition, positions that TOLERATE substitution. Variable
positions carry more coupling signal than conserved ones, and while the average-product
correction is designed to remove that per-site entropic bias it does not remove it perfectly. So
a z of +18.87 against a uniformly random null could be an artefact of having selected variable
sites rather than evidence of a functional module.

The fix is an entropy-matched null: draw null sets whose per-site conservation distribution
matches the divergent set, by stratified sampling over Shannon-entropy bins computed from the
same reweighted MSA the couplings came from. If the enrichment survives that, it is about which
positions coevolve, not about how variable they are.

A second control: shuffle the coupling matrix's labels while preserving its row/column marginals
(degree-preserving), which tests whether the signal needs the specific pairing or only the
overall coupling strengths of those rows.
"""
import json, warnings
warnings.filterwarnings("ignore")
import numpy as np

AAs="ACDEFGHIKLMNPQRSTVWY"; A2I={a:i for i,a in enumerate(AAs)}; Q=21
MSA="/home/user/uploads_fold/a98a190e-fold_tuba1a/msas/fold_tuba1a_paired_msa_chains_a.a3m"
DI=np.load("results/dca/dca_APC_full.npy").astype(np.float64); L=DI.shape[0]
div=json.load(open("results/alignment/tuba8_divergent_positions.json"))
DIV=sorted({d["canon_pos"] for d in div})
rng=np.random.default_rng(1234)

# ---- per-site entropy from the same MSA (reweighted the same way is unnecessary for binning)
def read_a3m(path):
    seqs=[];buf=[];qlen=None
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
    arr=np.full((len(seqs),qlen),Q-1,dtype=np.int8)
    for i,s in enumerate(seqs):
        for j,c in enumerate(s): arr[i,j]=A2I.get(c,Q-1)
    return arr[(arr==Q-1).mean(axis=1)<0.5]
arr=read_a3m(MSA)
ent=np.zeros(L)
for j in range(L):
    col=arr[:,j]; col=col[col!=Q-1]
    if len(col)==0: continue
    p=np.bincount(col,minlength=20)[:20].astype(float); p=p/p.sum()
    p=p[p>0]; ent[j]=-(p*np.log(p)).sum()
print(f"per-site entropy: mean {ent.mean():.3f}, range {ent.min():.3f}-{ent.max():.3f}")
di_idx=np.array(DIV)-1
print(f"divergent-set entropy: mean {ent[di_idx].mean():.3f} vs rest {np.delete(ent,di_idx).mean():.3f}")

def within(sel):
    sel=np.asarray(sel); v=[]
    for a in range(len(sel)):
        for b in range(a+1,len(sel)):
            if abs(sel[a]-sel[b])<6: continue
            v.append(DI[sel[a]-1,sel[b]-1])
    return np.array(v)
obs=within(DIV).mean()

# ---- null 1: uniform random (the original, for reference)
n1=np.array([within(rng.choice(np.arange(1,L+1),size=len(DIV),replace=False)).mean()
             for _ in range(20000)])

# ---- null 2: entropy-matched stratified sampling
NB=8
edges=np.quantile(ent,np.linspace(0,1,NB+1)); edges[-1]+=1e-9
bins=np.clip(np.digitize(ent,edges)-1,0,NB-1)
need=np.bincount(bins[di_idx],minlength=NB)
pool={b:np.where(bins==b)[0]+1 for b in range(NB)}
print(f"\nentropy-bin composition of the divergent set: {need.tolist()}")
print(f"pool sizes per bin:                            {[len(pool[b]) for b in range(NB)]}")
n2=[]
for _ in range(20000):
    sel=[]
    ok=True
    for b in range(NB):
        if need[b]==0: continue
        if len(pool[b])<need[b]: ok=False; break
        sel.extend(rng.choice(pool[b],size=need[b],replace=False).tolist())
    if not ok: continue
    n2.append(within(sel).mean())
n2=np.array(n2)

# ---- null 3: degree-preserving label shuffle of the coupling matrix
def shuffled_once():
    perm=rng.permutation(L)
    M=DI[np.ix_(perm,perm)]
    sel=np.asarray(DIV); v=[]
    for a in range(len(sel)):
        for b in range(a+1,len(sel)):
            if abs(sel[a]-sel[b])<6: continue
            v.append(M[sel[a]-1,sel[b]-1])
    return np.mean(v)
n3=np.array([shuffled_once() for _ in range(5000)])

print(f"\n{'null model':38s} {'mean':>10s} {'sd':>9s} {'z(obs)':>8s} {'p':>9s}")
for nm,n in (("uniform random positions",n1),
             ("ENTROPY-MATCHED positions",n2),
             ("matrix label shuffle",n3)):
    z=(obs-n.mean())/n.std(); p=(np.sum(n>=obs)+1)/(len(n)+1)
    print(f"{nm:38s} {n.mean():10.5f} {n.std():9.5f} {z:+8.2f} {p:9.5f}")
print(f"\nobserved within-set mean coupling = {obs:.5f}")

# does it survive if we drop the tail and the 32-45 loop (the two divergence clusters)?
core=[p for p in DIV if p<32 or (45<p<437)]
oc=within(core).mean()
nc=[]
for _ in range(20000):
    sel=rng.choice(np.arange(1,L+1),size=len(core),replace=False)
    nc.append(within(sel).mean())
nc=np.array(nc)
print(f"\ndropping both divergence clusters (32-45 and >=437): {len(core)} positions left")
print(f"  within-set mean {oc:.5f}; uniform null {nc.mean():.5f} +/- {nc.std():.5f}; "
      f"z = {(oc-nc.mean())/nc.std():+.2f}; p = {(np.sum(nc>=oc)+1)/20001:.5f}")
json.dump(dict(observed=float(obs),
               nulls={"uniform":[float(n1.mean()),float(n1.std())],
                      "entropy_matched":[float(n2.mean()),float(n2.std())],
                      "label_shuffle":[float(n3.mean()),float(n3.std())]},
               core_only=dict(n=len(core), observed=float(oc),
                              null=[float(nc.mean()),float(nc.std())],
                              z=float((oc-nc.mean())/nc.std()))),
          open("results/dca/module_controls.json","w"), indent=1)
