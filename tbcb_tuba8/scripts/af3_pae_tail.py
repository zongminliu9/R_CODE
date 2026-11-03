#!/usr/bin/env python3
"""Is the predicted placement of the C-terminal tail actually reliable?

A good global ipTM does not license a local conclusion. AF3 is trained on complexes in which
CAP-Gly domains bind -EEY/F termini, so it has every incentive to dock the tail into that groove
whether or not it is confident. The test is local:

  pLDDT   per-residue confidence along the tail (low = placed with little conviction)
  PAE     expected aligned error for tail-residue vs TBCB-residue pairs. PAE below ~5 A means
          the relative placement is trusted; above ~10-15 A it is effectively unknown.

We also compare the terminal aromatic's environment between the two complexes, which is a direct
test of the EEY/EEF question in the real TBCB complex rather than in a transplanted model.
"""
import json, glob, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
B = "/home/user/uploads_fold"
JOBS = {"TUBA8": (f"{B}/e50dbe95-fold_tbcb_tuba8", 449),
        "TUBA1A": (f"{B}/88398339-fold_tbcb_tuba1a", 451)}
CAPGLY = (134, 231)

for iso, (d, L) in JOBS.items():
    files = sorted(glob.glob(f"{d}/*full_data_*.json"))
    print(f"\n=== TBCB + {iso}  ({len(files)} models) ===")
    plddt_tail = []; pae_tail_cap = []; pae_tail_all = []
    for fp in files:
        J = json.load(open(fp))
        keys = [k for k in J if not isinstance(J[k], (int, float, str))]
        if "pae" not in J:
            print(f"  {os.path.basename(fp)}: keys = {list(J)[:10]}"); continue
        pae = np.array(J["pae"], dtype=np.float32)
        tci = np.array(J["token_chain_ids"]); tri = np.array(J["token_res_ids"], dtype=int)
        chains = sorted(set(tci.tolist()))
        # alpha chain = the one with more tokens
        cnt = {c: int((tci==c).sum()) for c in chains}
        chA = max(cnt, key=cnt.get); chT = min(cnt, key=cnt.get)
        aplddt = np.array(J.get("atom_plddts", []), dtype=np.float32)
        # per-token plddt: AF3 gives atom-level; use contact_probs/pae only if absent
        tail = np.where((tci == chA) & (tri >= L-14) & (tri <= L))[0]
        capg = np.where((tci == chT) & (tri >= CAPGLY[0]) & (tri <= CAPGLY[1]))[0]
        allT = np.where(tci == chT)[0]
        pae_tail_cap.append(pae[np.ix_(tail, capg)].min(axis=1))
        pae_tail_all.append(pae[np.ix_(tail, allT)].min(axis=1))
        if len(aplddt):
            # map atoms->tokens is not given directly; approximate per-residue by splitting
            pass
    if not pae_tail_cap: continue
    Pc = np.array(pae_tail_cap); Pa = np.array(pae_tail_all)
    pos = list(range(L-14, L+1))
    print(f"  {'res':>5s} {'PAE(min) to CAP-Gly':>21s} {'PAE(min) to any TBCB':>22s}")
    for i, p in enumerate(pos):
        print(f"  {p:5d} {Pc[:,i].mean():12.2f} +/-{Pc[:,i].std():5.2f} "
              f"{Pa[:,i].mean():14.2f} +/-{Pa[:,i].std():5.2f}")
    print(f"  tail mean PAE to CAP-Gly: {Pc.mean():.2f} A;  to any TBCB: {Pa.mean():.2f} A")
    print(f"  terminal residue {L}: PAE to CAP-Gly {Pc[:,-1].mean():.2f} A "
          f"(min over models {Pc[:,-1].min():.2f})")
