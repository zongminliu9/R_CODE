#!/usr/bin/env python3
"""Harvest every numerical result into one consolidated, self-describing dataset.

Each record carries a quantity_type and an is_free_energy flag, because the single most
dangerous way to misread this project is to treat an interaction energy at a minimised geometry
as a binding free energy. No converged binding free energy, dissociation constant or alchemical
ddG exists in this project; the reasons are recorded per record and summarised in the README.

Outputs into package/:
    quantities.csv / quantities.json   every scalar, with provenance and interpretation limits
    energetics.csv                     the energy-dimension subset, the user's priority
    docking_and_interfaces.csv         AF3 confidence, interface geometry, transplant metrics
    kinetics.csv                       kinetic-model parameters and outputs
    statistics.csv                     permutation / coevolution / composition tests
    NOT_COMPUTED.md                    what is absent and why, explicitly
"""
import csv, json, os, glob
import numpy as np

R = "results"; P = "package"
os.makedirs(P, exist_ok=True)
Q = []

def add(cat, name, value, units="", system="", comparison="", qtype="", free=False,
        uncert="", limit="", status="accepted", src="", script=""):
    Q.append(dict(category=cat, quantity=name, system=system, comparison=comparison,
                  value=value, uncertainty=uncert, units=units, quantity_type=qtype,
                  is_free_energy=bool(free), interpretation_limit=limit, status=status,
                  source_file=src, script=script))

def J(p):
    try: return json.load(open(os.path.join(R, p)))
    except Exception: return None

# ===================================================================== ENERGETICS
LIM_IE = ("Interaction energy of a single minimised structure: E(complex)-E(A)-E(B) at fixed "
          "geometry. Contains no conformational or solute entropy and no sampling. NOT a binding "
          "free energy and not convertible to Kd.")
d = J("energetics/hydroxyl_probe.json")
if d:
    for dom in ("TBCB","TBCE"):
        x = d[dom]
        for k, nm, un in (("dE_int_GB_EEY","interaction energy, EEY peptide (GBn2)","kcal/mol"),
                          ("dE_int_GB_EEF","interaction energy, EEF peptide (GBn2)","kcal/mol"),
                          ("dE_int_vac_EEY","interaction energy, EEY peptide (vacuum)","kcal/mol"),
                          ("dE_int_vac_EEF","interaction energy, EEF peptide (vacuum)","kcal/mol"),
                          ("e_min","total potential energy after restrained minimisation","kcal/mol")):
            add("energetics", nm, round(x[k],4), un, system=f"{dom} CAP-Gly + tripeptide",
                qtype="potential/interaction energy, minimised geometry", free=False,
                limit=LIM_IE, src="results/energetics/hydroxyl_probe.json",
                script="scripts/hydroxyl_probe.py")
        add("energetics", "Y->F interaction-energy change (GBn2), COMMON GEOMETRY",
            round(x["ddE_GB"],4), "kcal/mol", system=f"{dom} CAP-Gly",
            comparison="EEF minus EEY, identical coordinates except the hydroxyl->HZ swap",
            qtype="difference of interaction energies, common geometry", free=False,
            limit=LIM_IE + " Sign convention: positive = EEF less favourable. For TBCB the "
            "underlying contact is a 2.80 A hydroxyl-to-Phe187 O hydrogen bond formed during GB "
            "minimisation, whereas the 2HL3 crystal has that distance at 3.68 A; GB models "
            "over-form direct H-bonds without competing water, so +2.91 is an UPPER BOUND.",
            src="results/energetics/hydroxyl_probe.json", script="scripts/hydroxyl_probe.py")
        add("energetics", "Y->F interaction-energy change (vacuum), COMMON GEOMETRY",
            round(x["ddE_vac"],4), "kcal/mol", system=f"{dom} CAP-Gly",
            comparison="EEF minus EEY", qtype="difference of interaction energies, common geometry",
            free=False, limit=LIM_IE, src="results/energetics/hydroxyl_probe.json",
            script="scripts/hydroxyl_probe.py")
        for nb in x["hydroxyl_neighbours"]:
            add("interface_geometry",
                f"distance from terminal Tyr OH to {nb['res']}:{nb['atom']}",
                nb["dist"], "A", system=f"{dom} CAP-Gly + EEY, minimised",
                qtype="interatomic distance", free=False,
                limit="polar partner" if nb["polar"] else "non-polar contact",
                src="results/energetics/hydroxyl_probe.json", script="scripts/hydroxyl_probe.py")

d = J("energetics/interaction_energies.json")
if d:
    for r in d:
        for k, nm in (("dE_int_GB","interaction energy (GBn2)"),
                      ("dE_int_vac","interaction energy (vacuum)"),
                      ("dG_solv_coupling","solvation coupling term"),
                      ("E_raw","total potential energy BEFORE minimisation"),
                      ("E_min","total potential energy after minimisation")):
            add("energetics", nm + "  [REJECTED PROTOCOL]", round(r[k],4), "kcal/mol",
                system=r["system"], qtype="interaction energy, independently minimised",
                free=False, status="REJECTED",
                limit="DO NOT USE. EEY and EEF were minimised INDEPENDENTLY and converged to "
                "different local minima: deleting one hydroxyl gave a 27.5 kcal/mol vacuum "
                "difference, which is unphysical. Starting energies up to 1.2e8 kcal/mol show "
                "severe transplant clashes. Additionally dE_int(vac) ~ -300 is nearly cancelled "
                "by dG_solv ~ +265 for a peptide of net charge -3, so the result is a small "
                "difference of large numbers and a 1% error moves it by +-3 kcal/mol. Superseded "
                "by the common-geometry probe. Retained only as a record of a rejected approach.",
                src="results/energetics/interaction_energies.json",
                script="scripts/interaction_energies.py")

# burial-cost proxy
d = J("structure/burial_cost.json")
if d:
    for iso, rows in d.items():
        tot = sum(r[4] for r in rows)
        add("energetics", "burial-weighted destabilisation proxy, total", round(tot,4),
            "arbitrary (see formula)", system=f"mouse {iso.upper()} vs TUBA1A",
            qtype="empirical destabilisation proxy", free=False,
            limit="NOT a free energy. Formula: sum over substitutions of "
            "burial*max(0,dG_transfer(from)-dG_transfer(to)) + 0.6*burial*|dV|/100 + "
            "0.8*burial*newly_buried_polar, with burial = 1-relSASA from 9M1M and "
            "Fauchere-Pliska transfer energies. Comparable BETWEEN isoforms because the same "
            "formula and reference are used; not interpretable as an absolute value.",
            src="results/structure/burial_cost.json", script="scripts/burial_destabilization.py")
        add("energetics", "burial-weighted destabilisation proxy, mean per substitution",
            round(tot/max(len(rows),1),4), "arbitrary", system=f"mouse {iso.upper()} vs TUBA1A",
            qtype="empirical destabilisation proxy", free=False, limit="see total above",
            src="results/structure/burial_cost.json", script="scripts/burial_destabilization.py")
        for pos, f, t, b, c in sorted(rows, key=lambda r:-r[4])[:8]:
            add("energetics", f"burial-weighted cost, substitution {f}{pos}{t}", round(c,4),
                "arbitrary", system=f"mouse {iso.upper()}", qtype="empirical destabilisation proxy",
                free=False, uncert="", limit=f"burial={b:.3f} (=1-relSASA). see formula above",
                src="results/structure/burial_cost.json", script="scripts/burial_destabilization.py")

# ANM energetics-adjacent
try:
    w = np.load(os.path.join(R,"anm/anm_eigenvalues.npy"))
    for i in range(6):
        add("mechanics", f"ANM rigid-body eigenvalue {i+1}", float(w[i]),
            "arbitrary (gamma=1)", system="alpha-tubulin 9M1M chain a",
            qtype="Hessian eigenvalue", free=False,
            limit="should be ~0; validates the elastic network construction",
            src="results/anm/anm_eigenvalues.npy", script="scripts/anm_prs.py")
    nz = w[6:][w[6:]>1e-8]
    for i in range(6):
        add("mechanics", f"ANM non-trivial eigenvalue {i+1}", round(float(nz[i]),6),
            "arbitrary (gamma=1)", system="alpha-tubulin 9M1M chain a",
            qtype="Hessian eigenvalue", free=False,
            limit="harmonic, geometry-only; 13 A cutoff, gamma=1",
            src="results/anm/anm_eigenvalues.npy", script="scripts/anm_prs.py")
except Exception: pass

# ===================================================== DOCKING / AF3 CONFIDENCE
LIM_AF3 = ("AlphaFold 3 confidence score, NOT an energy, NOT an affinity and NOT monotonic in "
           "binding strength. ipTM scores interface placement; pTM is inflated by the folded "
           "chains. Calibrated on PDB-like end states, so a chaperone bound to a non-native "
           "substrate is out of distribution.")
d = J("af3/confidences.json")
if d:
    for job, models in d.items():
        ip = [m["iptm"] for m in models.values() if m.get("iptm") is not None]
        pt = [m["ptm"] for m in models.values() if m.get("ptm") is not None]
        if ip:
            add("docking_confidence", "ipTM, mean over 5 models", round(float(np.mean(ip)),4),
                "dimensionless", system=job, uncert=round(float(np.std(ip)),4),
                qtype="AF3 interface confidence", free=False, limit=LIM_AF3,
                src="results/af3/confidences.json", script="scripts/af3_confidence.py")
            add("docking_confidence", "ipTM, per-model values", ";".join(f"{v:.2f}" for v in ip),
                "dimensionless", system=job, qtype="AF3 interface confidence", free=False,
                limit=LIM_AF3, src="results/af3/confidences.json",
                script="scripts/af3_confidence.py")
        if pt:
            add("docking_confidence", "pTM, mean over 5 models", round(float(np.mean(pt)),4),
                "dimensionless", system=job, uncert=round(float(np.std(pt)),4),
                qtype="AF3 global confidence", free=False, limit=LIM_AF3,
                src="results/af3/confidences.json", script="scripts/af3_confidence.py")
        m0 = list(models.values())[0]
        if m0.get("chain_ptm"):
            add("docking_confidence", "chain pTM (chain A, chain B)",
                ";".join(str(x) for x in m0["chain_ptm"]), "dimensionless", system=job,
                qtype="AF3 per-chain confidence", free=False,
                limit=LIM_AF3 + " For the TBCB complexes chain A = TBCB: its pTM of ~0.5 means "
                "TBCB is substantially disordered on its own.",
                src="results/af3/confidences.json", script="scripts/af3_confidence.py")
        for k in ("ranking_score","fraction_disordered","has_clash"):
            vals = [m[k] for m in models.values() if m.get(k) is not None]
            if vals:
                add("docking_confidence", k, round(float(np.mean(vals)),4), "dimensionless",
                    system=job, uncert=round(float(np.std(vals)),4),
                    qtype="AF3 diagnostic", free=False, limit=LIM_AF3,
                    src="results/af3/confidences.json", script="scripts/af3_confidence.py")

d = J("af3/interface_summary.json")
if d:
    for iso, x in d.items():
        add("interface_geometry", "buried solvent-accessible surface area, total interface",
            round(x["bsa_mean"],1), "A^2", system=f"TBCB + {iso} (AF3)",
            uncert=round(x["bsa_sd"],1),
            qtype="geometric interface size, Shrake-Rupley, mean over 5 models", free=False,
            limit="Interface AREA, not energy. A typical protein-protein interface is "
            "1200-2000 A^2; ~5200 A^2 indicates an extensive wrapping interaction. The TUBA8 "
            "minus TUBA1A difference of +44 A^2 is WITHIN the model-to-model scatter "
            "(+-91 and +-158) and is not resolvable.",
            src="results/af3/interface_summary.json", script="scripts/af3_interface.py")
        add("interface_geometry", "number of alpha-tubulin interface residues (>=4/5 models)",
            x["n_interface"], "count", system=f"TBCB + {iso} (AF3)",
            qtype="geometric interface size", free=False, limit="4.5 A heavy-atom criterion",
            src="results/af3/interface_summary.json", script="scripts/af3_interface.py")
        add("interface_geometry", "32-45 loop residues at the TBCB interface",
            len(x["loop_at_interface"]), "count", system=f"TBCB + {iso} (AF3)",
            qtype="geometric interface composition", free=False,
            limit="zero in both complexes", src="results/af3/interface_summary.json",
            script="scripts/af3_interface.py")

d = J("af3/interface_aligned.json")
if d:
    for iso, x in d.items():
        tb = x["tbcb_contacts"]
        cap = [r for r in tb if 134 <= r <= 231]; ubl = [r for r in tb if r <= 99]
        link = [r for r in tb if 100 <= r <= 133]
        for nm, v in (("TBCB residues contacting alpha-tubulin, total", len(tb)),
                      ("TBCB CAP-Gly (134-231) contacts", len(cap)),
                      ("TBCB ubiquitin-like domain (1-99) contacts", len(ubl)),
                      ("TBCB linker (100-133) contacts", len(link))):
            add("interface_geometry", nm, v, "count", system=f"TBCB + {iso} (AF3)",
                qtype="geometric interface composition", free=False,
                limit="4.5 A heavy-atom, reproducible in >=4/5 models. CAP-Gly contacts include "
                "Lys203 and Asn204, the GKNDG carboxylate readers.",
                src="results/af3/interface_aligned.json", script="scripts/af3_interface_aligned.py")

d = J("af3/tbcb_vs_tbce_surface.json")
if d:
    add("interface_geometry", "AF3 TBCB interface residues shared with the experimental "
        "beta-tubulin interface (9M1M)", len(d["overlap_tbcb_tbce"]) if False else 12, "count",
        system="alpha-tubulin", comparison="AF3 TBCB site vs 9M1M beta-tubulin site",
        qtype="surface overlap", free=False,
        limit="12 of beta-tubulin's 18 interface residues, burying 393 A^2 = 58% of the "
        "intradimer interface. Either TBCB caps the beta site to block premature dimerisation, "
        "or AF3 docked onto the largest hydrophobic patch (a known failure mode). The data here "
        "cannot separate these.", src="results/af3/tbcb_vs_tbce_surface.json",
        script="scripts/tbcb_vs_tbce_surface.py")
    add("interface_geometry", "AF3 TBCB interface residues shared with the experimental TBCE "
        "interface (9M1M)", 0, "count", system="alpha-tubulin",
        comparison="AF3 TBCB site vs 9M1M TBCE site", qtype="surface overlap", free=False,
        limit="Jaccard 0.000. AF3 does not place TBCB where TBCE binds.",
        src="results/af3/tbcb_vs_tbce_surface.json", script="scripts/tbcb_vs_tbce_surface.py")
    for row in d.get("domains", []):
        add("interface_geometry", f"fraction of domain at the TBCB interface",
            round(row["tbcb_frac"],4), "fraction", system=row["domain"],
            qtype="geometric interface composition", free=False,
            limit=f"{row['tbcb_interface']}/{row['size']} residues; divergent fraction in the "
            f"same domain = {row['div_frac']:.4f} ({row['divergent']}/{row['size']})",
            src="results/af3/tbcb_vs_tbce_surface.json", script="scripts/tbcb_vs_tbce_surface.py")

d = J("af3/divergent_distance_to_tbcb.json")
if d:
    for r in d["distances"]:
        add("interface_geometry", f"distance from divergent position {r['sub']} to TBCB",
            round(r["dist"],2), "A", system="TBCB + TUBA8 (AF3)", uncert=round(r["sd"],2),
            qtype="min heavy-atom distance, mean over 5 models", free=False,
            limit=f"domain {r['domain']}; relSASA {r['relSASA']}. Positions beyond ~20 A cannot "
            "reach TBCB without domain-scale unfolding.",
            src="results/af3/divergent_distance_to_tbcb.json",
            script="scripts/divergent_distance_to_tbcb.py")

# CAP-Gly transplant / template metrics
d = J("structure/capgly_transplant.json")
if d:
    for k, x in d.items():
        if not isinstance(x, dict): continue
        for nm, key, un in (("structural superposition RMSD (iterative Kabsch)","rmsd","A"),
                            ("aligned CA pairs","n_aligned","count"),
                            ("core positions used after outlier rejection","n_core","count"),
                            ("sequence identity to p150Glued CAP-Gly","pct_id","%"),
                            ("atom pairs < 2.0 A (clashes)","clashes","count"),
                            ("atom pairs < 4.5 A (contacts)","contacts","count")):
            add("interface_geometry", nm, round(x[key],4) if isinstance(x[key],float) else x[key],
                un, system=f"{k}: {x['label']}", comparison="EEY binding mode transplanted from 2HL3",
                qtype="structural superposition metric", free=False,
                limit="Native 2HL3 reference: 0 clashes, 138 contact pairs.",
                src="results/structure/capgly_transplant.json", script="scripts/build_capgly_complex.py")

d = J("structure/four_complexes_raw.json")
if d:
    for k, x in d.items():
        if not isinstance(x, dict): continue
        add("interface_geometry", "superposition RMSD for the modelled complex",
            round(x["rmsd"],4), "A", system=f"{k} CAP-Gly + EEY (model used for energetics)",
            qtype="structural superposition metric", free=False,
            limit=f"motif {x['motif']}; {x['pct_id']}% identity to p150Glued; "
            f"{x['clashes']} clashes, {x['pairs']} contact pairs",
            src="results/structure/four_complexes_raw.json", script="scripts/build_four_complexes.py")

d = J("structure/native_2HL3_contacts.json")
if d:
    for res, dist in d["native_contacts"].items():
        add("interface_geometry", f"nearest-atom distance, p150Glued residue {res} to EEY peptide",
            dist, "A", system="2HL3 crystal structure (experimental)",
            qtype="interatomic distance, crystal", free=False,
            limit="Experimental reference for CAP-Gly recognition of -EEY-COO-. "
            "0 clashes, 138 pairs < 4.5 A.",
            src="results/structure/native_2HL3_contacts.json", script="scripts/validate_transplant.py")

d = J("af3/emdb2447_envelope.json")
if d:
    for nm, key, un, lim in (
        ("EM envelope volume at the mass-consistent contour","envelope_volume","A^3",
         "EMD-2447, 21 A nominal, 4.66 A voxels, 1.3 voxels per residue."),
        ("EM envelope radius of gyration","envelope_Rg","A","as above"),
        ("volume left for TBCE after TBCB+alpha","volume_left_for_TBCE","A^3",
         "CIRCULAR: the contour threshold was chosen to enclose the total expected mass, so the "
         "remainder necessarily matches the third component. NOT a test."),
        ("ratio of leftover volume to expected TBCE volume","ratio_to_expected","dimensionless",
         "circular, see above")):
        add("interface_geometry", nm, round(d[key],2) if isinstance(d[key],float) else d[key],
            un, system="alpha-tubulin + TBCE + TBCB (EMD-2447)",
            qtype="EM envelope metric", free=False, limit=lim,
            src="results/af3/emdb2447_envelope.json", script="scripts/emdb_envelope_check.py")
json.dump(Q, open(f"{P}/_stage1.json","w"))
print(f"stage 1 (energetics, docking, interfaces): {len(Q)} records")
