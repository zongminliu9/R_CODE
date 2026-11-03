# TBCB-dependent TUBA8 folding and isoform-selective stability

Computational investigation of why TBCB knockout reduces TUBA8 protein in cerebellum while
other α-tubulin isoforms, and other brain regions, are comparatively preserved.

## Layout

```
data/seq     49 reviewed UniProt sequences (mouse + human), FASTA
data/meta    UniProt records, sequence index, structure index  (provenance)
data/pdb     30 experimental structures, mmCIF        [gitignored, regenerable]
data/afdb    8 AlphaFold v6 models                    [gitignored, regenerable]
scripts      analysis and simulation code
results      numerical output (CSV + JSON)
report       written findings
```

Structural data is not committed because it is large and byte-for-byte regenerable:

```bash
python3 scripts/fetch_sequences.py
python3 scripts/fetch_structures.py
```

## Scripts

| script | what it does |
|---|---|
| `fetch_sequences.py` | resolve + download reviewed UniProt entries for the tubulin/cofactor system |
| `fetch_structures.py` | download the working structure set from RCSB and AlphaFold DB |
| `survey_structures.py` / `find_capgly.py` | locate relevant PDB entries by search, with metadata |
| `inspect_capgly_templates.py` | report chain/ligand composition of CAP-Gly complexes |
| `align_isoforms.py` | exact Needleman–Wunsch alignment; divergent-position enumeration |
| `map_divergence_structure.py` | identify chains by sequence; set up structural mapping |
| `divergence_vs_interfaces.py` | per-residue SASA in 3 contexts, interface burial, cofactor distances |
| `enrichment_test.py` | permutation tests for structural-class enrichment; density scan |
| `maturation_model.py` | multistate maturation ODE model with shared limiting cofactor pools |
| `selectivity_analysis.py` | null model (Q0) and one-mechanism-at-a-time decomposition (Q1) |
| `regionality_sensitivity.py` | demand scan and demand × TBCE grid (Q2) |
| `global_sensitivity.py` | broad-prior sensitivity and phenotype-region characterisation (Q3) |
| `fast_steady.py` | Newton steady-state solver with integration fallback |
| `bench_openmm.py`, `bench_real.py`, `bench_tune.py` | measured MD throughput on this hardware |

## Findings

See `report/01_findings_to_date.md`. Headline results: TUBA8 is a ~2–20× outlier in family
divergence; its 45 substitutions are statistically *depleted* at every functional interface and
concentrated in two disordered, exposed regions (the H1–S2/K40 loop, 12 substitutions in
residues 32–45, p = 5×10⁻⁵; and the C-terminal tail); and in the kinetic model the only
parameter able to reproduce the observed selective loss is the rate of TBCB-independent capture
by TBCE.

## Hardware note

This workspace ran on 4 CPU cores with no GPU. Measured OpenMM throughput is in
`report/01_findings_to_date.md` §4. Converged protein–protein binding free energies, PMFs and
alchemical ΔΔG values are **not** obtainable at that throughput and are not reported; the
scripts and prepared inputs are structured so they can be run unchanged on GPU hardware.
