# What is NOT in this package, and why

Read this before the data files. The brief asked for "all binding-energy and affinity numbers
with no omissions". The honest answer is that **this project contains no binding free energies,
no dissociation constants and no alchemical ΔΔG values**, and the energy-dimension numbers it
does contain are not substitutes for them. Treating them as substitutes is the one mistake that
would invalidate a report built on this package.

---

## 1. Absent quantities

| requested / expected | present? | why not |
|---|---|---|
| ΔG_bind(TBCB, TUBA8 / TUBA1A / TUBA1B / TUBA4A) | **no** | needs a converged protein–protein PMF or double-decoupling: 30–50 windows × 20–100 ns = 1–5 μs per system. Measured throughput on the available hardware (no GPU) makes this decades. |
| K_d, or ΔG° = RT ln(K_d/C°) | **no** | follows only from a converged standard binding free energy, which does not exist here. Nothing in `quantities.csv` may be converted to a K_d. |
| ΔΔG_bind between isoforms | **no** | same reason. |
| alchemical ΔΔG_fold (FEP / TI) for any mutation | **no** | ~20 λ windows × 10 ns × 2 legs ≈ 400 ns per mutation. Inputs are built and validated in `gpu/`; the calculation has not been run. |
| umbrella-sampling or metadynamics PMFs | **no** | infeasible at the measured throughput. |
| MM/PBSA | **no** | only MM/GBSA was attempted, and it was rejected (§3). |
| explicit-solvent MD trajectories | **no** | α-tubulin monomer runs at 3.9 ns/day on this machine; the originally specified 4 isoforms × 3 replicas × 20 ns needs 62 days with an implicit-solvent shortcut already applied. |
| docking scores (Vina, HADDOCK, ZDOCK, Rosetta…) | **no** | **no docking program was run at any point.** The complexes come from AlphaFold 3 co-folding, which is structure prediction, not docking, and its ipTM is a confidence score with no energy dimension. Do not report any number in this package as a "docking score". |
| per-residue free-energy decomposition | **no** | requires converged ensembles. |
| transition-state barriers ΔG‡ from simulation | **no** | no reactive or enhanced sampling was performed. Barrier-like statements in the kinetic model are parameters, not computed barriers. |
| TBCB↔TBCE transfer free-energy surface | **no** | no structure of a TBCB–TBCE–α-tubulin complex at usable resolution exists (EMD-2447 is 21 Å, ~1.3 voxels per residue), and the sampling is infeasible here. |
| CCT/TRiC or TBCD/TBCC/ARL2 energetics | **no** | not attempted. |

## 2. What the energy-dimension numbers in this package actually are

Every record in `quantities.csv` carries `quantity_type` and `is_free_energy`. **`is_free_energy`
is `False` for all 512 records.** The energy-dimension records fall into three kinds:

**(a) Interaction energies at a minimised geometry** — `E(complex) − E(A) − E(B)` evaluated at one
fixed structure. No conformational entropy, no solute entropy, no sampling, no ensemble. Useful
only for comparing chemically similar systems that share a binding mode, which is why they were
used for the EEY-versus-EEF question and nothing else. Units kcal/mol, which makes them look
like free energies; they are not.

**(b) An empirical destabilisation proxy** — the burial-weighted cost used to compare isoforms.
Its units are arbitrary and the formula is stated in the record. It is comparable *between*
isoforms because the formula and reference are identical; an absolute value means nothing.

**(c) Elastic-network eigenvalues and perturbation-response couplings** — harmonic, geometric,
with γ = 1, so in arbitrary units. Not energies of any physical process.

## 3. One calculation was performed and then rejected

20 records are flagged `status = REJECTED`, from `scripts/interaction_energies.py`. They are kept
so that nobody re-derives them or, worse, uses them. Three independent reasons:

1. EEY and EEF were minimised **independently** and converged to different local minima: deleting
   a single hydroxyl produced a 27.5 kcal/mol vacuum difference, which is unphysical.
2. Starting energies reached 1.2 × 10⁸ kcal/mol, i.e. severe clashes in the transplanted models.
3. For a peptide of net charge −3, ΔE_int(vacuum) ≈ −300 kcal/mol is nearly cancelled by
   ΔG_solv ≈ +265, so the result is a small difference of large numbers and a 1% error in either
   term moves it by ±3 kcal/mol — larger than the effect being measured.

They are superseded by the common-geometry probe (`scripts/hydroxyl_probe.py`), whose records are
`status = accepted` but still **not** free energies.

## 4. Numbers that exist but have no discriminating power

Reported because their own built-in controls invalidate them, and an analyst who does not know
that might build on them:

- **AlphaFold pLDDT of the 32–45 loop.** TUBA1A, TUBA1B and TUBA4A have an *identical* 32–45
  sequence yet span 66.09–76.81 (range 10.72), while TUBA8 deviates from their mean by 0.76.
- **Relative SASA of the 32–45 loop.** The same three identical-sequence isoforms span
  0.379–0.526; TUBA8 is 0.456.
- **The EMD-2447 volume bookkeeping.** The contour threshold was chosen to enclose the total
  expected mass, so "the leftover volume matches TBCE" follows by construction. It is circular and
  is not a test.
- **Buried-SASA difference between the two AF3 complexes.** +44 Å² for TUBA8, against
  model-to-model scatter of ±91 and ±158 Å². Not resolvable.

## 5. Kinetic-model parameters are placeholders

The 14 baseline rate constants are physically plausible values, not measurements: no experimental
k_on or k_off for TBCB–α-tubulin exists. **Every conclusion drawn from the model comes from
parameter *ratios* and from the 50 000-draw sensitivity sweep spanning ±1 decade on every
baseline and ±1.5 decades on every isoform ratio — not from the baseline values themselves.**
Quoting a baseline as a measured rate would be a misreading.

## 6. Two structural caveats that limit several conclusions

- **9M1M contains no TBCB.** Every enrichment and burial statistic computed on 9M1M constrains
  the TBCD / TBCE / TBCC / ARL2 surfaces of the *late* maturation complex, not the TBCB interface.
- **The AF3 complexes model a fully folded substrate** (chain pTM 0.91 for α-tubulin), whereas
  TBCB acts on a quasi-native intermediate released from CCT. And the predicted TBCB interface
  covers 12 of the 18 α-tubulin residues that β-tubulin uses in the mature heterodimer, burying
  58% of the intradimer interface — which is either the real capping mechanism or a known AF3
  failure mode, and the data here cannot separate them. TBCB's own chain pTM is 0.51.
