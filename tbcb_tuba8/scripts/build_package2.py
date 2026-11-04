#!/usr/bin/env python3
"""Stage 2: kinetics, statistical tests, sequence composition, coevolution, benchmarks.
Then write every output file into package/."""
import csv, json, os
import numpy as np
R = "results"; P = "package"
Q = json.load(open(f"{P}/_stage1.json"))

def add(cat, name, value, units="", system="", comparison="", qtype="", free=False,
        uncert="", limit="", status="accepted", src="", script=""):
    Q.append(dict(category=cat, quantity=name, system=system, comparison=comparison,
                  value=value, uncertainty=uncert, units=units, quantity_type=qtype,
                  is_free_energy=bool(free), interpretation_limit=limit, status=status,
                  source_file=src, script=script))
def J(p):
    try: return json.load(open(os.path.join(R, p)))
    except Exception: return None

# ======================================================================= KINETICS
BASE = dict(k_fold=1e-2, k_degU=1e-3, kon_B=1.0, koff_B=1e-2, k_degF=1e-3, kon_Edir=1e-2,
            k_tr=1.0, k_degB=1e-5, k_degE=1e-4, k_asm=1.0, k_pol=1e-1, k_depol=1e-2,
            k_degH=1e-4, k_degM=1e-5)
UNITS = {"kon_B":"1/(uM*s)","kon_Edir":"1/(uM*s)","k_tr":"1/(uM*s)"}
LIM_K = ("Reference parameter value of the multistate maturation model. These are PHYSICALLY "
         "PLAUSIBLE PLACEHOLDERS, not measured rate constants: no experimental kon/koff for "
         "TBCB-alpha-tubulin exists. Conclusions are drawn from parameter RATIOS and from the "
         "50000-draw sensitivity analysis spanning +-1 decade on every baseline, not from these "
         "values.")
for k, v in BASE.items():
    add("kinetics", f"model baseline parameter {k}", v, UNITS.get(k,"1/s"),
        system="multistate maturation model", qtype="rate constant (reference value)",
        free=False, limit=LIM_K, src="scripts/maturation_model.py",
        script="scripts/maturation_model.py")

d = J("kinetics/Q0_null_model.json")
if d:
    for r in d:
        add("kinetics", "selectivity (log2 FC_TUBA8 - log2 FC_TUBA1A), NULL MODEL",
            round(r["sel"],6), "log2 units",
            system=f"s_TUBA8/s_TUBA1A = {r['ratio']}",
            comparison="TBCB knockout vs wild type",
            qtype="model output", free=False,
            limit="NULL MODEL: both isoforms given IDENTICAL rate constants, differing only in "
            "synthesis flux. Selectivity is exactly 0.000000 at every ratio down to 0.01, so "
            "mass-action competition for a shared limiting chaperone pool produces NO isoform "
            "selectivity. Combined with the proteomics result that TUBA8 is quantified in all 60 "
            "TMT channels across all five brain regions, an expression-floor explanation of the "
            "regional phenotype is excluded and a TUBA8-specific rate difference is required.",
            src="results/kinetics/Q0_null_model.json", script="scripts/selectivity_analysis.py")
        add("kinetics", "fold change FC_TUBA8 (KO/WT), NULL MODEL", round(r["fc_TUBA8"],6),
            "ratio", system=f"s_TUBA8/s_TUBA1A = {r['ratio']}", qtype="model output",
            free=False, limit="see selectivity record",
            src="results/kinetics/Q0_null_model.json", script="scripts/selectivity_analysis.py")

d = J("kinetics/Q1_mechanism_decomposition.json")
if d:
    for mech, rows in d.items():
        for r in rows:
            add("kinetics", "selectivity, ONE MECHANISM AT A TIME", round(r["sel"],4),
                "log2 units", system=mech,
                comparison=f"TUBA8/TUBA1A parameter ratio = {r['ratio']}",
                qtype="model output", free=False,
                limit="Single parameter varied around one baseline, so interaction effects are "
                "invisible here; see the Q3 global analysis, which revises the reading of "
                "k_degF. Strongest reachable selectivity per mechanism: kon_Edir -3.30, "
                "kon_B -0.46, k_degF -0.32, k_fold 0.000, k_tr +0.03 (wrong sign), k_asm 0.000.",
                src="results/kinetics/Q1_mechanism_decomposition.json",
                script="scripts/selectivity_analysis.py")
            add("kinetics", "FC_TUBA8 (KO/WT), one mechanism at a time", round(r["fc_TUBA8"],4),
                "ratio", system=mech, comparison=f"ratio = {r['ratio']}",
                qtype="model output", free=False, limit="see selectivity record",
                src="results/kinetics/Q1_mechanism_decomposition.json",
                script="scripts/selectivity_analysis.py")

d = J("kinetics/Q2_regionality_demand.json")
if d:
    for r in d:
        add("kinetics", "selectivity vs synthesis demand", round(r["sel"],4), "log2 units",
            system=f"demand multiplier = {r['demand']}",
            comparison="TUBA8 carries ONE difference: kon_Edir reduced 30x",
            qtype="model output", free=False,
            limit="Selectivity peaks near demand ~ 1 and vanishes at high demand, the opposite "
            "of a naive 'cerebellum has the highest tubulin demand' expectation. The 2-D scan "
            "shows the |sel| ridge tracking demand ~ TBCE_tot.",
            src="results/kinetics/Q2_regionality_demand.json",
            script="scripts/regionality_sensitivity.py")
d = J("kinetics/Q2_grid_demand_TBCE.json")
if d:
    for r in d:
        add("kinetics", "selectivity, demand x TBCE grid", round(r["sel"],4), "log2 units",
            system=f"demand={r['demand']}, TBCE_tot={r['TBCE']}", qtype="model output",
            free=False, limit="ridge follows demand ~ TBCE_tot",
            src="results/kinetics/Q2_grid_demand_TBCE.json",
            script="scripts/regionality_sensitivity.py")

d = J("kinetics/Q3_global_sensitivity.json")
if d:
    add("kinetics", "usable draws in the global sensitivity sweep", d["n_usable"], "count",
        system="broad log-uniform prior", qtype="sampling diagnostic", free=False,
        limit=f"of {d['N']} drawn; phenotype reproduced in {d['n_phenotype']} "
        f"({100*d['n_phenotype']/d['n_usable']:.2f}%). Prior: every baseline rate +-1 decade, "
        "every TUBA8/TUBA1A ratio +-1.5 decades, demand/TBCB_tot/TBCE_tot over 3 decades.",
        src="results/kinetics/Q3_global_sensitivity.json",
        script="scripts/global_sensitivity_fast.py")
    for r in d["spearman"]:
        add("kinetics", f"Spearman rho of selectivity vs {r['param']}", round(r["rho"],4),
            "dimensionless", system="broad prior, 49728 draws", uncert=f"p={r['p']:.3e}",
            qtype="global sensitivity coefficient", free=False,
            limit="lr_ = log10 TUBA8/TUBA1A ratio; lb_ = log10 baseline value; log_ = context "
            "variable. Top three: kon_Edir +0.445, k_degF -0.293, kon_B -0.255. k_fold is "
            "irrelevant everywhere in the prior (rho -0.007, p=0.11).",
            src="results/kinetics/Q3_global_sensitivity.json",
            script="scripts/global_sensitivity_fast.py")
    for r in d["phenotype_region"]:
        add("kinetics", f"phenotype-region median shift in {r['param']}", round(r["shift"],4),
            "log10 units", system="draws with sel<-1.5, FC8<0.5, FC1a>0.75",
            uncert=f"p={r['p']:.3e}", qtype="phenotype-region characterisation", free=False,
            limit="Median in the phenotype region minus median in the rest. The phenotype needs "
            "a CONJUNCTION: k_degF +0.879, kon_Edir -1.014, kon_B +0.689, koff_B -0.227. "
            "11.9% of phenotype draws need no kon_Edir deficit at all, so the observation does "
            "not uniquely imply one.",
            src="results/kinetics/Q3_global_sensitivity.json",
            script="scripts/global_sensitivity_fast.py")

# ===================================================================== STATISTICS
d = J("alignment/pairwise_identity_mouse.json")
if d:
    for pair, x in d.items():
        add("sequence", "pairwise sequence identity", x["pct_id"], "%", system=pair.replace("|"," vs "),
            qtype="exact Needleman-Wunsch global alignment, BLOSUM62, affine gaps -11/-1",
            free=False, limit=f"{x['n_substitutions']} substitutions, {x['n_gap_cols']} gap "
            f"columns, {x['identities']}/{x['aligned_cols']} identities",
            src="results/alignment/pairwise_identity_mouse.json", script="scripts/align_isoforms.py")

d = J("structure/enrichment_permutation.json")
if d:
    for r in d:
        add("statistics", f"divergent-position enrichment at: {r['cls']}", r["observed"], "count",
            system="9M1M late cofactor complex", uncert=f"expected {r['expected']}",
            comparison=f"ratio {r['ratio']}",
            qtype="permutation test, 1e5 draws over 437 resolved residues", free=False,
            limit=f"p_enrich={r['p_enrich']:.5f}, p_deplete={r['p_deplete']:.5f}. "
            "NOTE: 9M1M contains NO TBCB, so this constrains the TBCD/TBCE/TBCC/ARL2 surfaces "
            "of the late complex, not the TBCB interface.",
            src="results/structure/enrichment_permutation.json", script="scripts/enrichment_test.py")

d = J("dca/dca_validation.json")
if d:
    for r in d:
        add("statistics", f"DCA contact precision, top {r['n']} pairs, |i-j|>={r['sep']}",
            round(r["precision"],4), "fraction", system="alpha-tubulin, 9M1M Cb-Cb < 8 A",
            uncert=f"baseline {r['baseline']:.4f}", comparison=f"enrichment {r['enrichment']:.2f}x",
            qtype="method validation", free=False,
            limit="Validates that the coevolutionary couplings carry real information before any "
            "allosteric claim rests on them. 60.7% on the top L/4 long-range pairs is within the "
            "normal working range for mean-field DCA.",
            src="results/dca/dca_validation.json", script="scripts/dca2.py")
d = J("dca/dca_meta_full.json")
if d:
    for k, nm in (("M","sequences used"),("Meff","effective sequences after 80% reweighting"),
                  ("Meff_over_L","Meff/L")):
        add("statistics", f"DCA {nm}", round(d[k],3) if isinstance(d[k],float) else d[k], "count",
            system="AF3 paired MSA of TUBA1A", qtype="sampling diagnostic", free=False,
            limit="Meff/L = 2.42 is below the ideal >5 for mean-field DCA, so couplings are "
            "noisier than optimal; the contact-precision validation is what justifies their use.",
            src="results/dca/dca_meta_full.json", script="scripts/dca2.py")
d = J("dca/allostery.json")
if d:
    for r in d["loop_coupling"]:
        add("statistics", f"32-45 loop non-contacting coevolutionary coupling to: {r['target']}",
            round(r["loop"],6), "APC-corrected coupling", system="alpha-tubulin",
            uncert=f"null mean {r['null_mean']:.6f}", comparison=f"percentile {r['percentile']:.1f}",
            qtype="coevolutionary coupling, Cb-Cb > 12 A and |i-j| >= 6", free=False,
            limit=f"p={r['p']:.4f}, null = all other contiguous 14-residue windows "
            f"(n={r['n_null_windows']}), {r['n_pairs']} qualifying pairs. The loop is at the "
            "0.0 percentile for the TBCB interface, i.e. LESS coupled than any other window.",
            src="results/dca/allostery.json", script="scripts/dca_allostery.py")
    m = d["module"]
    add("statistics", "within-set coevolutionary coupling of the 45 divergent positions",
        round(m["mean"],6), "APC-corrected coupling", system="alpha-tubulin",
        uncert=f"uniform null {m['null_mean']:.6f} +- {m['null_sd']:.6f}",
        qtype="coevolutionary module test", free=False,
        limit=f"p={m['p']:.5f} against a uniform null, but see the entropy-matched control: "
        "divergent positions are more variable (entropy 1.425 vs 0.694) and that absorbs most of "
        "the naive z. Entropy-matched z = +5.76; dropping both spatial clusters leaves 26 "
        "scattered substitutions at z = +4.44.",
        src="results/dca/allostery.json", script="scripts/dca_allostery.py")
d = J("dca/module_controls.json")
if d:
    for nm, v in d["nulls"].items():
        z = (d["observed"]-v[0])/v[1]
        add("statistics", f"coevolving-module z-score against the {nm} null", round(z,3),
            "z", system="45 TUBA8-divergent positions",
            uncert=f"null {v[0]:.5f} +- {v[1]:.5f}", qtype="controlled permutation test",
            free=False, limit="The entropy-matched null is the one that matters; the uniform "
            "null is inflated by the divergent set being more variable.",
            src="results/dca/module_controls.json", script="scripts/dca_module_control.py")
    c = d["core_only"]
    add("statistics", "coevolving-module z-score, scattered substitutions only",
        round(c["z"],3), "z", system=f"{c['n']} divergent positions outside both clusters",
        uncert=f"null {c['null'][0]:.5f} +- {c['null'][1]:.5f}",
        qtype="controlled permutation test", free=False,
        limit="Both the 32-45 loop and the C-terminal tail removed, so this is the scattered "
        "surface substitutions alone. They still coevolve.",
        src="results/dca/module_controls.json", script="scripts/dca_module_control.py")

d = J("anm/loop32_45_coupling.json")
if d:
    for r in d:
        add("mechanics", f"32-45 loop ANM mechanical coupling to: {r['target']}",
            round(r["coupling"],8), "PRS units", system="alpha-tubulin 9M1M chain a",
            uncert=f"null mean {r['null_mean']:.8f}", comparison=f"percentile {r['percentile']:.1f}",
            qtype="elastic-network perturbation response", free=False,
            limit=f"p={r['p']:.4f}. HARMONIC and GEOMETRY-ONLY: blind to entropy-driven and "
            "side-chain-mediated coupling, so a null result here does not exclude allostery. The "
            "DCA test is the one that does. Note also that with a symmetric covariance the PRS "
            "matrix is symmetric, so effectiveness and sensitivity are identical and carry no "
            "directional information.",
            src="results/anm/loop32_45_coupling.json", script="scripts/anm_coupling_test.py")

d = J("structure/plddt_loop32_45.json")
if d:
    for iso, x in d.items():
        add("statistics", "AlphaFold pLDDT, 32-45 loop mean", round(x["loop"],2), "pLDDT",
            system=iso, uncert=f"rest of chain {x['rest']:.2f}",
            qtype="model confidence", free=False,
            limit="NO DISCRIMINATING POWER. The three isoforms with IDENTICAL 32-45 sequence span "
            "66.09-76.81 (range 10.72) while TUBA8 deviates from their mean by only 0.76. The "
            "loop is low-confidence in every alpha-tubulin.",
            src="results/structure/plddt_loop32_45.json", script="scripts/plddt_loop.py")
d = J("structure/loop32_45_properties.json")
if d:
    for iso, x in d.items():
        if "loop_rel_sasa" in x:
            add("statistics", "relative SASA of the 32-45 loop", x["loop_rel_sasa"], "fraction",
                system=iso, qtype="geometric", free=False,
                limit="NO DISCRIMINATING POWER. Three isoforms with identical loop sequence span "
                "0.379-0.526, a spread exceeding any TUBA8 deviation.",
                src="results/structure/loop32_45_properties.json",
                script="scripts/loop32_45_properties.py")

d = J("structure/folding_units.json")
if d:
    rco = np.array(d["rel_contact_order"]); res = d["residues"]
    add("statistics", "chain-average relative contact order", round(float(rco.mean()),4),
        "dimensionless", system="alpha-tubulin 9M1M chain a",
        qtype="topological folding-rate predictor", free=False,
        limit="0.0653 is LOW for a 437-residue globular protein (typically 0.08-0.15), so the "
        "topology is an easy fold. Yet tubulin is obligately chaperone-dependent, which means its "
        "chaperone requirement is not a topological folding-difficulty problem -- matching the "
        "kinetic model's independent finding that k_fold has no effect on selectivity.",
        src="results/structure/folding_units.json", script="scripts/folding_units.py")
    r2i = {r:i for i,r in enumerate(res)}
    loop = [r2i[p] for p in range(32,46) if p in r2i]
    add("statistics", "relative contact order, 32-45 loop", round(float(rco[loop].mean()),4),
        "dimensionless", system="alpha-tubulin", uncert=f"rest {np.delete(rco,loop).mean():.4f}",
        qtype="topological", free=False,
        limit="p(low)=0.0000. The LOWEST contact order in the protein, and the loop has ZERO "
        "contacts to any other folding unit (p=0.0000): it forms early and autonomously and does "
        "not participate in the late cooperative steps where a chaperone is required.",
        src="results/structure/folding_units.json", script="scripts/folding_units.py")

# composition
add("sequence", "substitutions to Ser or Thr, TUBA1A->TUBA8", 12, "count",
    system="45 divergent positions (43 scorable)", uncert="expected 4.9",
    qtype="binomial test against background composition", free=False,
    limit="background S+T frequency in TUBA1A = 11.3%, binomial p = 0.0023. To Ser alone: 7 "
    "observed vs 2.0 expected, p = 0.0035. Net +6 hydroxyl side chains. NOT spatially clustered "
    "(mean pairwise CA 29.54 A vs null 27.96+-2.44, p=0.74) and NOT a distinguishable "
    "coevolutionary sub-module (z=+0.97, p=0.17).",
    src="scripts/ser_substitutions.py", script="scripts/ser_substitutions.py")
add("sequence", "net charge change contributed by the 45 divergent positions", 0, "e",
    system="TUBA1A -> TUBA8", qtype="formal charge count", free=False,
    limit="Whole protein: TUBA1A -24, TUBA8 -23. Exactly zero net change across 45 "
    "substitutions is a strong constraint signature.",
    src="scripts/ser_substitutions.py", script="scripts/ser_substitutions.py")
add("sequence", "burial-cost placement test: observed vs position-randomised", -3.35, "z",
    system="TUBA8 substitutions", uncert="null 22.27 +- 1.50, observed 17.23",
    qtype="permutation test", free=False,
    limit="p(higher)=0.9991. Keeping TUBA8's substitution identities but randomising their "
    "positions RAISES the burial cost, so the real placement avoids buried positions. Third "
    "buffering layer alongside zero net charge change and low contact order (p=0.0036).",
    src="results/structure/burial_cost.json", script="scripts/burial_destabilization.py")

# ===================================================================== BENCHMARKS
for nm, atoms, cond, nsday in (
    ("CAP-Gly + EEY peptide (2HL3)", 1227, "GBn2 implicit, 2.0 nm cutoff, 4 fs", 6.2),
    ("CAP-Gly + EEY peptide (2HL3)", 1227, "GBn2, 1.0 nm cutoff, 4 fs", 20.3),
    ("CAP-Gly + EEY peptide (2HL3)", 1227, "GBn2, 1.2 nm, HMR, 5 fs", 17.1),
    ("CAP-Gly + EEY peptide (2HL3)", 1227, "HCT implicit, 1.2 nm, HMR, 5 fs, 4 threads", 26.0),
    ("CAP-Gly + EEY peptide (2HL3)", 1227, "HCT, 1 thread (x4 parallel = 35.0 aggregate)", 8.75),
    ("mouse TBCB CAP-Gly (1WHG)", 1650, "GBn2, 2.0 nm, 4 fs", 3.6),
    ("alpha-tubulin monomer (1JFF)", 6375, "GBn2, 1.0 nm, HMR, 5 fs", 3.88)):
    add("benchmark", "OpenMM throughput", nsday, "ns/day", system=f"{nm}, {atoms} atoms",
        comparison=cond, qtype="measured simulation throughput", free=False,
        limit="4 x Intel Xeon 2.10 GHz, 15 GB RAM, NO GPU, OpenMM 8.6.1 CPU platform. "
        "Consequence: Simulation A1 as originally specified (4 isoforms x 3 replicas x 20 ns of "
        "alpha-tubulin monomer) needs 62 days even with the implicit-solvent shortcut; the "
        "specified explicit ff19SB/OPC is a further 10-30x slower. Protein-protein PMF and "
        "per-mutation FEP are 10^3-10^4 GPU-days.",
        src="report/01_findings_to_date.md", script="scripts/bench_real.py, scripts/bench_tune.py")

os.makedirs(P, exist_ok=True)
cols = ["category","quantity","system","comparison","value","uncertainty","units",
        "quantity_type","is_free_energy","interpretation_limit","status","source_file","script"]
with open(f"{P}/quantities.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in Q: w.writerow({k: r.get(k,"") for k in cols})
json.dump(Q, open(f"{P}/quantities.json","w"), indent=1)
for cat, fn in (("energetics","energetics.csv"),
                ("docking_confidence","docking_and_interfaces.csv"),
                ("interface_geometry","docking_and_interfaces.csv"),
                ("kinetics","kinetics.csv"),
                ("statistics","statistics.csv"),
                ("mechanics","statistics.csv"),
                ("sequence","statistics.csv"),
                ("benchmark","benchmarks.csv")):
    rows = [r for r in Q if r["category"] == cat]
    mode = "a" if os.path.exists(f"{P}/{fn}") and fn in ("docking_and_interfaces.csv","statistics.csv") else "w"
    with open(f"{P}/{fn}", mode, newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        if mode == "w": w.writeheader()
        for r in rows: w.writerow({k: r.get(k,"") for k in cols})
os.remove(f"{P}/_stage1.json")
print(f"total records: {len(Q)}")
from collections import Counter
for c, n in sorted(Counter(r["category"] for r in Q).items(), key=lambda x:-x[1]):
    print(f"  {c:22s} {n:5d}")
nfe = sum(1 for r in Q if r["is_free_energy"])
print(f"\nrecords flagged as free energies: {nfe}")
print(f"records flagged REJECTED:        {sum(1 for r in Q if r['status']=='REJECTED')}")
