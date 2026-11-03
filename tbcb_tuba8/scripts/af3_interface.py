#!/usr/bin/env python3
"""The TBCB-alpha-tubulin interface, from the AF3 predictions.

This is the interface that has no experimental structure and that every earlier conclusion in
this project had to leave unmapped. The question it settles is whether TUBA8's divergent
positions -- and in particular the 32-45 loop -- actually touch TBCB.

Analysed over all 5 models per complex so that a contact is only called if it is reproducible,
not if it appears in one seed. Interface burial is computed as dSASA by the same Shrake-Rupley
protocol used for 9M1M, so the numbers are comparable with the earlier tables.
"""
import json, glob, os, csv, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}
MAXSASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
           'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}
B = "/home/user/uploads_fold"
JOBS = {"TUBA8": f"{B}/e50dbe95-fold_tbcb_tuba8", "TUBA1A": f"{B}/88398339-fold_tbcb_tuba1a"}
LOOP = set(range(32, 46))
div = json.load(open("results/alignment/tuba8_divergent_positions.json"))
DIV = {d["canon_pos"] for d in div}

def chains_of(cif):
    arr = pdbx.get_structure(pdbx.CIFFile.read(cif), model=1, use_author_fields=True)
    arr = arr[struc.filter_amino_acids(arr)]
    info = {}
    for ch in np.unique(arr.chain_id):
        sub = arr[arr.chain_id == ch]
        ca = sub[sub.atom_name == "CA"]
        info[ch] = (ca.array_length(), "".join(AA3.get(n,"X") for n in ca.res_name[np.argsort(ca.res_id)]))
    return arr, info

def res_sasa(ctx, chain):
    s = struc.sasa(ctx, vdw_radii="Single", point_number=200)
    out = {}
    for ch, rid, v in zip(ctx.chain_id, ctx.res_id, s):
        if ch != chain or np.isnan(v): continue
        out[int(rid)] = out.get(int(rid), 0.0) + float(v)
    return out

summary = {}
for iso, d in JOBS.items():
    models = sorted(glob.glob(f"{d}/*_model_*.cif"))
    print(f"\n=== TBCB + {iso}  ({len(models)} models) ===")
    arr0, info = chains_of(models[0])
    for ch,(n,sq) in info.items():
        print(f"    chain {ch}: {n} res   {sq[:34]}...{sq[-10:]}")
    # TBCB is the shorter chain
    chT = min(info, key=lambda c: info[c][0]); chA = max(info, key=lambda c: info[c][0])
    print(f"    -> TBCB = chain {chT} ({info[chT][0]} aa), alpha-tubulin = chain {chA} ({info[chA][0]} aa)")

    per_model = []
    for m in models:
        arr = pdbx.get_structure(pdbx.CIFFile.read(m), model=1, use_author_fields=True)
        arr = arr[struc.filter_amino_acids(arr)]
        A = arr[arr.chain_id == chA]; T = arr[arr.chain_id == chT]
        sa_free = res_sasa(A, chA); sa_cplx = res_sasa(A + T, chA)
        st_free = res_sasa(T, chT); st_cplx = res_sasa(A + T, chT)
        dA = {r: sa_free[r]-sa_cplx.get(r,0.0) for r in sa_free}
        dT = {r: st_free[r]-st_cplx.get(r,0.0) for r in st_free}
        bsa = sum(dA.values()) + sum(dT.values())
        dmin = np.linalg.norm(A.coord[:,None,:]-T.coord[None,:,:], axis=-1)
        contactsA = {int(A.res_id[i]) for i in np.where(dmin.min(axis=1) < 4.5)[0]}
        per_model.append(dict(model=os.path.basename(m), bsa=bsa, dA=dA, dT=dT,
                              contactsA=contactsA))
        print(f"    {os.path.basename(m):34s} buried SASA {bsa:7.1f} A^2  "
              f"alpha interface residues {len(contactsA)}")

    # reproducible interface: residue counted if it contacts in >=4 of 5 models
    cnt = {}
    for pm in per_model:
        for r in pm["contactsA"]: cnt[r] = cnt.get(r,0)+1
    robust = sorted([r for r,c in cnt.items() if c >= 4])
    dA_mean = {}
    for pm in per_model:
        for r,v in pm["dA"].items(): dA_mean.setdefault(r,[]).append(v)
    dA_mean = {r: float(np.mean(v)) for r,v in dA_mean.items()}
    aa1 = {int(r): AA3.get(n,"X") for r,n in zip(arr0[arr0.chain_id==chA].res_id,
                                                 arr0[arr0.chain_id==chA].res_name)}
    print(f"\n    reproducible alpha-tubulin interface ({len(robust)} residues, >=4/5 models):")
    print("      " + " ".join(f"{aa1.get(r,'X')}{r}{'*' if r in DIV else ''}" for r in robust))
    inloop = [r for r in robust if r in LOOP]
    indiv  = [r for r in robust if r in DIV]
    print(f"    32-45 loop residues at the interface: {inloop if inloop else 'NONE'}")
    print(f"    TUBA8-divergent positions at the interface: {indiv if indiv else 'NONE'}")
    print(f"    mean buried SASA over models: {np.mean([p['bsa'] for p in per_model]):.1f} "
          f"+/- {np.std([p['bsa'] for p in per_model]):.1f} A^2")
    summary[iso] = dict(chain_tbcb=chT, chain_alpha=chA,
                        bsa_mean=float(np.mean([p['bsa'] for p in per_model])),
                        bsa_sd=float(np.std([p['bsa'] for p in per_model])),
                        robust_interface=robust, loop_at_interface=inloop,
                        divergent_at_interface=indiv,
                        dSASA_alpha=dA_mean,
                        n_interface=len(robust))
json.dump({k:{kk:vv for kk,vv in v.items() if kk!="dSASA_alpha"} for k,v in summary.items()},
          open("results/af3/interface_summary.json","w"), indent=1)
json.dump({k:v["dSASA_alpha"] for k,v in summary.items()},
          open("results/af3/interface_dSASA.json","w"), indent=1)

print("\n=== TUBA8 vs TUBA1A interface comparison ===")
s8, s1 = summary["TUBA8"], summary["TUBA1A"]
print(f"  buried SASA:  TUBA8 {s8['bsa_mean']:7.1f} +/- {s8['bsa_sd']:.1f}   "
      f"TUBA1A {s1['bsa_mean']:7.1f} +/- {s1['bsa_sd']:.1f}   "
      f"diff {s8['bsa_mean']-s1['bsa_mean']:+7.1f} A^2")
print(f"  interface residues: TUBA8 {s8['n_interface']}   TUBA1A {s1['n_interface']}")
a, b = set(s8["robust_interface"]), set(s1["robust_interface"])
print(f"  shared {len(a&b)}; TUBA8-only {sorted(a-b)}; TUBA1A-only {sorted(b-a)}")
