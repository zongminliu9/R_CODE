#!/usr/bin/env python3
"""Three candidate loci for TUBA8's TBCB dependence have now been excluded: the functional
interfaces (depleted, p=0.004), allosteric coupling (n.s.), and C-terminal tail recognition
(refuted by TUBA4A, which cannot be read by CAP-Gly at all yet is preserved). What remains is
the H1-S2/K40 loop, residues 32-45, where 12 of 45 substitutions sit (p=5e-5).

The kinetic model says the viable mechanism is NOT slower folding (k_fold has zero effect) but
faster clearance of the chaperone-free quasi-native state (k_degF, the strongest discriminator
of the phenotype region). A loop that packs poorly would do exactly that. So we ask whether
TUBA8's loop is predicted to be less ordered and more hydrophobically exposed than TUBA1A's.

Cheap, model-based evidence only: AlphaFold pLDDT (per-residue confidence, which tracks local
disorder) and solvent exposure of hydrophobic atoms in the deposited models. No sampling.
"""
import json, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V'}
MAXSASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
           'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,
      'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
HYDROPHOBIC = set("AVILMFWCY")

MODELS = {"TUBA8":"Q9JJZ2","TUBA1A":"P68369","TUBA1B":"P05213","TUBA4A":"P68368"}
LOOP = list(range(32, 46))

idx = json.load(open("data/meta/sequence_index.json"))
GENE = {"TUBA8":"mouse_Tuba8","TUBA1A":"mouse_Tuba1a","TUBA1B":"mouse_Tuba1b","TUBA4A":"mouse_Tuba4a"}

print("=== sequence character of the 32-45 segment ===")
print(f"{'isoform':9s} {'32-45 sequence':>16s} {'Gly':>4s} {'Pro':>4s} {'chg':>4s} "
      f"{'<KD hydropathy>':>16s}")
seqrows = {}
for iso, acc in MODELS.items():
    s = idx[GENE[iso]]["sequence"]
    seg = s[31:45]
    ng = seg.count("G"); npro = seg.count("P")
    chg = sum(1 for c in seg if c in "KR") - sum(1 for c in seg if c in "DE")
    kd = np.mean([KD[c] for c in seg])
    print(f"{iso:9s} {seg:>16s} {ng:4d} {npro:4d} {chg:+4d} {kd:16.2f}")
    seqrows[iso] = dict(segment=seg, n_gly=ng, n_pro=npro, net_charge=chg, kd=float(kd))

print("\n=== AlphaFold pLDDT: 32-45 loop versus the rest of the chain ===")
print(f"{'isoform':9s} {'pLDDT loop':>11s} {'pLDDT rest':>11s} {'delta':>7s} "
      f"{'min in loop':>12s} {'rank of loop mean':>18s}")
out = {}
for iso, acc in MODELS.items():
    f = pdbio.PDBFile.read(f"data/afdb/{acc}.pdb")
    arr = f.get_structure(model=1)
    ca = arr[arr.atom_name == "CA"]
    o = np.argsort(ca.res_id); ca = ca[o]
    pl = f.get_b_factor(model=1) if hasattr(f, "get_b_factor") else None
    # biotite exposes B-factors via the 'b_factor' annotation when present
    bf = arr.b_factor if hasattr(arr, "b_factor") else None
    if bf is None:
        print(f"{iso}: no B-factor column"); continue
    per = {}
    for rid, b in zip(arr.res_id, bf):
        per.setdefault(int(rid), []).append(float(b))
    plddt = {r: float(np.mean(v)) for r, v in per.items()}
    loopv = [plddt[r] for r in LOOP if r in plddt]
    restv = [v for r, v in plddt.items() if r not in LOOP]
    # where does the loop mean rank among all 14-residue windows?
    L = max(plddt); keys = sorted(plddt)
    wins = [np.mean([plddt[k] for k in keys[i:i+14]]) for i in range(len(keys)-13)]
    rank = int(np.sum(np.array(wins) < np.mean(loopv)))
    print(f"{iso:9s} {np.mean(loopv):11.2f} {np.mean(restv):11.2f} "
          f"{np.mean(loopv)-np.mean(restv):+7.2f} {min(loopv):12.2f} "
          f"{f'{rank}/{len(wins)}':>18s}")
    out[iso] = dict(plddt_loop=float(np.mean(loopv)), plddt_rest=float(np.mean(restv)),
                    plddt_loop_min=float(min(loopv)), window_rank=rank, n_windows=len(wins))

print("\n=== exposed hydrophobic surface of the 32-45 loop (AlphaFold models) ===")
print(f"{'isoform':9s} {'loop SASA':>10s} {'hydrophobic':>12s} {'frac hydroph':>13s} "
      f"{'rel.SASA loop':>14s}")
for iso, acc in MODELS.items():
    f = pdbio.PDBFile.read(f"data/afdb/{acc}.pdb")
    arr = f.get_structure(model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    sasa = struc.sasa(arr, vdw_radii="Single", point_number=300)
    tot = hyd = 0.0; relnum = relden = 0.0
    aa1 = {int(r): AA3.get(n, "X") for r, n in zip(arr.res_id, arr.res_name)}
    for rid, nm, el, s in zip(arr.res_id, arr.atom_name, arr.element, sasa):
        if int(rid) not in LOOP or np.isnan(s): continue
        tot += float(s)
        if el == "C" or (aa1.get(int(rid)) in HYDROPHOBIC and el in ("C","S")):
            hyd += float(s)
    for r in LOOP:
        if r in aa1: relden += MAXSASA.get(aa1[r], 200)
    print(f"{iso:9s} {tot:10.1f} {hyd:12.1f} {hyd/max(tot,1e-9):13.3f} "
          f"{tot/max(relden,1e-9):14.3f}")
    out.setdefault(iso, {}).update(loop_sasa=round(tot,1), loop_hydrophobic_sasa=round(hyd,1),
                                   loop_frac_hydrophobic=round(hyd/max(tot,1e-9),3),
                                   loop_rel_sasa=round(tot/max(relden,1e-9),3))
for iso in out: out[iso].update(seqrows.get(iso, {}))
json.dump(out, open("results/structure/loop32_45_properties.json","w"), indent=1)
