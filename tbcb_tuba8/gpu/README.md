# GPU package — the two calculations this project actually needs

Everything in the main report that could be settled on 4 CPU cores has been settled. What is
left are two free-energy questions, and they are left because they need sampling, not because
they need more analysis. This directory is set up to run them unchanged on Sherlock or Marlowe.

## What is being computed, and why these two

The surviving mechanism (report 03, §5) is **folded-state stability of the chaperone-free
α-tubulin monomer** — the `k_degF` term that the 50 000-draw sensitivity analysis identified as
the strongest discriminator of the phenotype region. Two candidates remain:

**A. I114S** — the one substitution in TUBA8 that is not buffered. 91% buried, losing
1.84 kcal/mol of hydrophobic burial and 78 Å³ of side-chain volume, 2.4× the next-largest insult.
A single-atom-scale change is exactly the regime where alchemical free energy converges, so this
is the calculation with the best return.

**B. the 32–45 loop swap** — strictly conserved in all five other mouse α-tubulins and replaced
only in TUBA8, a perfect correlation with the phenotype. A 14-residue segment substitution is not
tractable alchemically, so it is addressed by replica-exchange sampling of local unfolding
instead of by a single ΔΔG.

A third calculation (TBCB binding free energy) is **deliberately not prioritised**: the AF3
complexes show the two interfaces are superimposable, with none of TUBA8's 38 core substitutions
in contact, so there is no structural hypothesis for a ΔΔG_bind to test. The tail-only relative
calculation is included as a cheap check because it is small.

## Jobs

| id | calculation | method | system size | estimate (1×A100) |
|---|---|---|---|---|
| A1 | ΔΔG_fold(I114S), folded leg | alchemical FEP, 20 λ, MBAR | ≈70 k atoms | ≈40 GPU-h |
| A2 | ΔΔG_fold(I114S), unfolded reference | FEP on Ace-GIS-Nme tripeptide | ≈4 k atoms | ≈3 GPU-h |
| A3 | controls: I110V, G232S, A294S, A270V | same protocol | ≈70 k atoms | ≈40 GPU-h each |
| B1 | TUBA8 vs TUBA8[loop→1A] local unfolding | REST2, 16 replicas | ≈70 k atoms | ≈200 GPU-h |
| B2 | same for TUBA1A vs TUBA1A[loop→8] | REST2, 16 replicas | ≈70 k atoms | ≈200 GPU-h |
| C1 | tail Y→F relative binding FE | FEP in complex and free | ≈120 k / 20 k atoms | ≈60 GPU-h |

A1+A2+A3 is the minimum useful run (≈170 GPU-h) and answers whether I114S destabilises TUBA8
enough to matter. B is the larger commitment.

## Running

```bash
module load python/3.12 cuda/12.4            # Sherlock; on Marlowe see slurm_marlowe.sbatch
conda env create -f environment.yml && conda activate tbcb
python 00_build_systems.py --all             # once, writes systems/
sbatch slurm_sherlock.sbatch A1              # or: A2 A3 B1 B2 C1
python 90_analyse_fep.py --job A1            # MBAR, overlap diagnostics, bootstrap CI
```

`00_build_systems.py` needs `../data/` and `../results/af3/` from this repository, and the AF3
model files; it writes nothing outside `systems/`.

## Protocol choices, and why

- **ff19SB with OPC water**, 310 K, 150 mM NaCl, as the original plan specified. ff19SB's
  CMAP-corrected backbone matters here because the question is local unfolding propensity.
- **Hydrogen mass repartitioning to 3 amu with a 4 fs timestep**, which is validated for ff19SB
  and roughly doubles throughput.
- **20 λ windows** for the Ile→Ser transformation, with electrostatics and sterics decoupled in
  separate stages and a soft-core for the sterics, because a buried hydrophobic-to-polar change
  has a large vdW component that single-stage schedules converge badly.
- **MBAR with overlap diagnostics**: report the overlap matrix and refuse to quote a ΔΔG if
  adjacent-window overlap drops below 0.03 anywhere. Hysteresis is checked by running the
  transformation in both directions.
- **Three independent replicas per leg from different velocity seeds**, with the ΔΔG uncertainty
  taken from the spread across replicas rather than from MBAR's analytical error, which is
  known to be optimistic.
- **REST2 rather than plain MD for the loop**, because the quantity of interest is a rare local
  unfolding event and plain MD at 310 K will not sample it in any affordable time.

## What would count as an answer

For A: a converged ΔΔG_fold(I114S) with a replica-spread uncertainty. If TUBA8 is destabilised by
≳1.5 kcal/mol relative to TUBA1A at this position, the `k_degF` mechanism has a concrete physical
origin and the kinetic model's required +0.879 dex in clearance rate becomes quantitatively
checkable (a 1.5 kcal/mol stability loss corresponds to roughly a 12-fold shift in the
folded/unfolded ratio at 310 K, which is the right order).

For B: a difference in the population of locally-unfolded loop states between the native and
swapped sequences, which is the microscopic version of the same claim.

A null result in both would mean the phenotype is not explained by monomer stability either, and
the remaining possibilities would be post-translational (the loop carries K40, whose acetylation
TUBA8 cannot undergo) or cell-type-specific rather than biophysical.
