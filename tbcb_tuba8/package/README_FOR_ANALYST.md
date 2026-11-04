# TBCB-dependent TUBA8 stability — complete computational package

Everything produced in this investigation, as self-describing tables plus the scripts that made
them. **Read `NOT_COMPUTED.md` first.** It states what is absent — there are no binding free
energies, no K_d values, no alchemical ΔΔG and no docking scores — and what the energy-dimension
numbers actually are.

## The question

TBCB knockout markedly reduces TUBA8 protein in cerebellum while other α-tubulin isoforms, and
TUBA8 in hippocampus, cortex, striatum and spinal cord, are comparatively preserved. Proteomics
confirms TUBA8 is quantified in all 60 TMT channels across all five regions, so the regional
restriction is not an expression floor.

## Files

| file | contents |
|---|---|
| `NOT_COMPUTED.md` | **read first**: absent quantities, interpretation limits, rejected results |
| `quantities.csv` / `.json` | all 512 scalar results, with `quantity_type`, `is_free_energy`, `interpretation_limit`, `status`, and the script that produced each |
| `energetics.csv` | the 73 energy-dimension records (20 of them `REJECTED`) |
| `docking_and_interfaces.csv` | 134 records: AF3 confidence, interface geometry, superposition metrics, distances |
| `kinetics.csv` | 221 records: model parameters and all four sweeps |
| `statistics.csv` | 77 records: permutation tests, coevolution, composition, mechanics |
| `benchmarks.csv` | 7 measured throughputs, which is why the free-energy calculations are absent |
| `MANIFEST.md` | data provenance: every sequence, structure and map, with accession and source |

In the repository alongside this package:

| path | contents |
|---|---|
| `../report/00_MASTER_SUMMARY.md` | the mechanism argument end to end, with every number |
| `../report/01…04_*.md` | the four working reports in the order the work was done, including the conclusions that were later overturned |
| `../scripts/` | 44 scripts; every number traces to one |
| `../results/` | raw outputs, including the full DCA coupling matrix and the 50 000-draw sweep |
| `../gpu/` | validated, unrun inputs for the two remaining free-energy calculations |
| `../data/` | sequences and provenance; structures are gitignored and regenerable by the fetch scripts |

## Reproducing

```bash
python3 scripts/fetch_sequences.py      # 49 reviewed UniProt entries
python3 scripts/fetch_structures.py     # 30 PDB entries + 8 AlphaFold models
```

Then any analysis script. Each is self-contained and writes to `results/`. Python 3.13 with
numpy, scipy, biotite, biopython, openmm, pdbfixer, mdtraj, MDAnalysis, networkx, matplotlib.

## How to read the conclusions

The work proceeded by elimination, and several early conclusions were overturned by later
evidence. The four working reports are kept in chronological order **including the parts that
turned out to be wrong**, because the corrections are part of the evidence:

- report 01 concluded the C-terminal tail could not be a TBCB contact because it is disordered in
  every experimental structure. Report 03 overturned this: the AF3 complexes put the tail at the
  interface with the terminal aromatic docked at PAE 3.6–3.8 Å.
- report 01 and 02 treated the elastic-network null result as excluding long-range allostery.
  Report 03 corrected the overstatement and replaced it with a coevolutionary test that is
  genuinely capable of detecting such coupling, and does not.
- report 02 read the burial proxy as showing TUBA8 uniquely destabilised. Report 03's placement
  permutation showed the opposite: the substitutions are buffered.

`00_MASTER_SUMMARY.md` carries the final state of each claim.

## Package contents as shipped

```
README_FOR_ANALYST.md       this file
NOT_COMPUTED.md             READ FIRST — absent quantities and interpretation limits
MANIFEST.md                 data provenance: every sequence, structure, map and AF3 job
quantities.csv / .json      all 512 scalar results, self-describing
energetics.csv              73 energy-dimension records (20 flagged REJECTED)
docking_and_interfaces.csv  134 records: AF3 confidence, interfaces, superpositions, distances
kinetics.csv                221 records: model parameters and all four sweeps
statistics.csv              77 records: permutation, coevolution, composition, mechanics
benchmarks.csv              7 measured throughputs
reports/                    00_MASTER_SUMMARY.md plus the four chronological working reports
raw_results/                every raw output, including the full 451x451 DCA coupling matrices,
                            the 437x437 PRS matrix and the 50000-draw sensitivity sweep (.npz)
scripts/                    all 44 analysis scripts; every number traces to one
gpu/                        validated, unrun inputs for the two remaining free-energy jobs
models/                     the CAP-Gly + peptide complexes built and minimised here
sequences/                  49 reviewed UniProt FASTA files plus the index
```

Experimental structures (30 mmCIF) and AlphaFold models (8 PDB) are **not** included because
they are large and byte-for-byte regenerable by `scripts/fetch_structures.py`. The AF3
predictions were supplied by the user and are not redistributed here.

## Column meanings in `quantities.csv`

| column | meaning |
|---|---|
| `category` | energetics / docking_confidence / interface_geometry / kinetics / statistics / mechanics / sequence / benchmark |
| `quantity` | what was measured |
| `system` | which molecular system or parameter setting |
| `comparison` | what it is being compared against, where applicable |
| `value`, `uncertainty`, `units` | the number; uncertainty may be an sd, a p-value or a null-model summary |
| `quantity_type` | **the physical nature of the number** — read this before using it |
| `is_free_energy` | `False` for every record in this project |
| `interpretation_limit` | what the number may and may not be used for |
| `status` | `accepted` or `REJECTED` (20 records) |
| `source_file`, `script` | provenance |
