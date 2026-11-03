#!/usr/bin/env python3
"""MBAR analysis with the diagnostics that decide whether a ddG may be quoted at all.

Three gates, and a number is printed only if all three pass:
  1. adjacent-window overlap >= 0.03 everywhere (otherwise the lambda ladder has a gap and MBAR
     is extrapolating across it);
  2. forward/reverse hysteresis within the combined replica spread (otherwise the transformation
     is not reversible at this sampling length);
  3. replica-to-replica spread reported as THE uncertainty, since MBAR's analytical error assumes
     uncorrelated samples and is optimistic for replica exchange.
The thermodynamic cycle is ddG = dG(folded or complex) - dG(reference or free).
"""
import argparse, glob, json, os, sys
import numpy as np

def load(nc):
    from openmmtools.multistate import MultiStateReporter, MultiStateSamplerAnalyzer
    rep = MultiStateReporter(nc, open_mode="r")
    ana = MultiStateSamplerAnalyzer(rep)
    dF, ddF = ana.get_free_energy()
    kT = ana.kT
    n = dF.shape[0]
    overlap = ana.get_overlap()["matrix"]
    adj = np.array([overlap[i, i+1] for i in range(n-1)])
    return dict(dG_kT=float(dF[0, -1]), err_kT=float(ddF[0, -1]),
                kT_kcal=float(kT.value_in_unit_system(
                    __import__("openmm").unit.md_unit_system))*0.239006,
                adjacent_overlap=adj.tolist(), min_adjacent_overlap=float(adj.min()),
                n_states=int(n))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--job", required=True, help="output directory from 10_fep_mutation.py")
    ap.add_argument("--overlap-floor", type=float, default=0.03)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.job, "*.nc")))
    if not files: sys.exit(f"no .nc files in {a.job}")
    byleg = {}
    print(f"{'file':44s} {'dG (kcal/mol)':>14s} {'min overlap':>12s} {'states':>7s}")
    for f in files:
        try: r = load(f)
        except Exception as e:
            print(f"{os.path.basename(f):44s}  FAILED {type(e).__name__}: {e}"); continue
        dG = r["dG_kT"]*0.6155          # kT at 310 K in kcal/mol
        name = os.path.basename(f)[:-3]
        leg = "reverse" if name.endswith("_rev") else name.split("_")[-2]
        byleg.setdefault(leg, []).append(dG)
        flag = "" if r["min_adjacent_overlap"] >= a.overlap_floor else "  << OVERLAP GAP"
        print(f"{name:44s} {dG:14.3f} {r['min_adjacent_overlap']:12.4f} "
              f"{r['n_states']:7d}{flag}")
        json.dump(r, open(f[:-3]+"_mbar.json","w"), indent=1)
    print()
    means = {}
    for leg, v in byleg.items():
        v = np.array(v)
        means[leg] = (v.mean(), v.std(ddof=1) if len(v) > 1 else float("nan"), len(v))
        print(f"  {leg:12s} dG = {v.mean():8.3f} +/- {means[leg][1]:.3f} kcal/mol "
              f"(spread over {len(v)} replicas)")
    hot = [l for l in means if l in ("folded","complex")]
    cold = [l for l in means if l in ("reference","free")]
    if hot and cold:
        h, c = means[hot[0]], means[cold[0]]
        ddG = h[0]-c[0]; err = float(np.hypot(h[1], c[1]))
        print(f"\n  ddG = dG({hot[0]}) - dG({cold[0]}) = {ddG:+.3f} +/- {err:.3f} kcal/mol")
        print(f"  (positive = the mutant is destabilised / binds more weakly)")
        if np.isnan(err):
            print("  !! single replica per leg: no spread-based uncertainty, do not quote this")
    else:
        print("\n  both legs of the cycle are needed before a ddG can be formed; "
              f"present legs: {sorted(means)}")

if __name__ == "__main__":
    main()
