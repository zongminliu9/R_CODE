# Local execution protocol — intermediate-state physics for TUBA8 / TBCB

Hand-off for local sessions. The cloud container could not run production dynamics (four
restarts; no trajectory survived), so **nothing in this protocol has been executed**. Everything
below is written to be run as-is.

**Scientific target.** TBCB does not act on native α-tubulin. It acts on a quasi-native
intermediate released from CCT. Every structural conclusion in `report/` therefore describes the
*end state*, and the open question is what the **intermediate ensemble** looks like, whether
TUBA8's differs from TUBA1A's, and whether TBCB binds those intermediates differently. That is
what Phases 2–4 are for. Phases 5–6 are the rigorous free energies.

**Read before interpreting anything**: `package/NOT_COMPUTED.md`.

---

## Phase 0 — environment

```bash
conda create -n tubulin -c conda-forge -y \
  python=3.11 gromacs=2024 openmm=8.2 openmmtools=0.23 pymbar=4 perses \
  ambertools=23 plumed=2.9 py-plumed gmx_mmpbsa mdtraj mdanalysis \
  deeptime scikit-learn pdbfixer biotite biopython parmed numpy scipy matplotlib
conda activate tubulin
```

```bash
gmx --version | head -3
python -c "import openmm; print(openmm.Platform.getPlatformByName('CUDA').getName())"
python -c "import perses, openmmtools, pymbar, deeptime; print('ok')"
```

Separate env for docking, which conflicts with the above:

```bash
conda create -n dock -c conda-forge -y python=3.10 haddock3 autodock-vina openbabel pdb-tools
```

Resource planning. One A100 or H100 gives roughly:

| system | atoms (ff19SB/OPC, 1.2 nm pad) | ns/day on A100 |
|---|---|---|
| α-tubulin monomer | ≈ 68 000 | 150–250 |
| TBCB alone | ≈ 42 000 | 300–400 |
| TBCB + α-tubulin | ≈ 118 000 | 90–150 |
| CAP-Gly + tripeptide | ≈ 26 000 | 450–600 |

Every wall-time figure below assumes 150 ns/day for the monomer and 110 for the complex. Scale
by your measured throughput — **measure it first** with a 1 ns test, do not trust these.

---

## Phase 1 — systems

Six α-tubulin monomers and three complexes. The mutants are what make the comparison
interpretable: without them a TUBA8-versus-TUBA1A difference cannot be attributed to any
particular substitution.

| id | construct | why |
|---|---|---|
| `a1a` | TUBA1A wild type | reference |
| `a8` | TUBA8 wild type | the phenotype isoform |
| `a1a_I114S` | TUBA1A with I114S | isolates the one unbuffered substitution (91% buried, 1.84 kcal/mol hydrophobicity and 78 Å³ lost, 2.4× the next-largest insult) |
| `a8_S114I` | TUBA8 with S114I | reciprocal; if it rescues, I114S is the determinant |
| `a8_loop1A` | TUBA8 with 32–45 → `PDGQMPSDKTIGGG` | the loop is strictly conserved in all five other isoforms and replaced only in TUBA8 |
| `a1a_loop8` | TUBA1A with 32–45 → `ADGTFGTQASKIND` | reciprocal |
| `tbcb` | mouse TBCB (Q9D1E6) | its AF3 chain pTM is 0.51, so sample it alone before trusting any complex |
| `tbcb_a8` | TBCB + TUBA8 (AF3 model 0) | |
| `tbcb_a1a` | TBCB + TUBA1A (AF3 model 0) | |

Starting coordinates: AlphaFold models `data/afdb/P68369.pdb` (TUBA1A) and `Q9JJZ2.pdb` (TUBA8);
complexes from the user's AF3 jobs. Build mutants with PyMOL or `pdb4amber` + `tleap`, then
**minimise before anything else** — the AF3 complex as deposited has an RMS atomic force of
**106.5 kcal/mol/Å** and a potential energy of **+31 810 kcal/mol** (measured, in
`results/decomposition/`), because AF3 does not optimise against a force field.

```bash
# per system, from a cleaned PDB
gmx pdb2gmx -f ${ID}.pdb -o ${ID}_proc.gro -p ${ID}.top -ff amber19sb -water tip4p -ignh
gmx editconf -f ${ID}_proc.gro -o ${ID}_box.gro -bt dodecahedron -d 1.2
gmx solvate  -cp ${ID}_box.gro -cs tip4p.gro -o ${ID}_solv.gro -p ${ID}.top
gmx grompp -f mdp/ions.mdp -c ${ID}_solv.gro -p ${ID}.top -o ions.tpr -maxwarn 2
echo SOL | gmx genion -s ions.tpr -o ${ID}_ions.gro -p ${ID}.top \
  -pname NA -nname CL -neutral -conc 0.15
gmx grompp -f mdp/em.mdp   -c ${ID}_ions.gro -p ${ID}.top -o em.tpr
gmx mdrun -deffnm em
gmx grompp -f mdp/nvt.mdp  -c em.gro   -r em.gro   -p ${ID}.top -o nvt.tpr   # 500 ps, posre
gmx mdrun -deffnm nvt -nb gpu -pme gpu -update gpu
gmx grompp -f mdp/npt.mdp  -c nvt.gro  -r nvt.gro  -p ${ID}.top -o npt.tpr   # 2 ns, posre off ramp
gmx mdrun -deffnm npt -nb gpu -pme gpu -update gpu
```

`mdp` settings that matter: `integrator = md`, `dt = 0.002` (or 0.004 with
`define = -DFLEXIBLE` off and hydrogen mass repartitioning via `gmx pdb2gmx -heavyh`),
`constraints = h-bonds`, `cutoff-scheme = Verlet`, `rvdw = rcoulomb = 1.0`, `coulombtype = PME`,
`tcoupl = v-rescale` with `ref-t = 310`, `pcoupl = C-rescale` at 1 bar,
`DispCorr = EnerPres`. **α-tubulin carries GTP·Mg²⁺ at the N-site** — keep it, parameterise with
the Meagher/Carlson triphosphate set in `amber19sb`, and verify the Mg coordination survives
equilibration. Dropping the nucleotide changes the N-domain's stability and invalidates the
whole intermediate analysis.

Sanity gates before proceeding: backbone RMSD plateau < 3 Å over the last 1 ns of NPT, no
drifting Mg²⁺, total charge zero, and density 1.00 ± 0.01 g/cm³.

---

## Phase 2 — generate the intermediate ensemble (the core of the brief)

Three routes. Run **(a)** always; **(b)** gives the free-energy profile; **(c)** is the cheap way
to seed Phase 3.

### (a) REST2 — unbiased population of near-native through partially unfolded

Solute tempering scales only the protein's terms, so water stays at 310 K and the replica count
scales with solute size rather than box size.

```bash
# 24 replicas, effective solute temperature 310 -> 520 K, geometric ladder
python scripts_local/make_rest2_ladder.py --n 24 --tmin 310 --tmax 520 --top ${ID}.top
mpirun -np 24 gmx_mpi mdrun -deffnm rest2_${ID} -multidir rep{00..23} \
  -replex 400 -hrex -plumed plumed_rest2.dat -nb gpu -pme gpu
```

200 ns per replica = 4.8 µs aggregate per system. On 4×A100, ≈ 4 days per system. Run `a8`,
`a1a`, `a8_loop1A` and `a1a_loop8` — that is the comparison the loop hypothesis needs.

Collective variables to record for every frame (PLUMED):

```
# plumed_rest2.dat
WHOLEMOLECULES ENTITY0=1-6900
q_all:   CONTACTMAP ... SWITCH={RATIONAL R_0=0.45 NN=8 MM=12} SUM   # global native contacts
q_nterm: CONTACTMAP ... # residues 1-205 only
q_loop:  CONTACTMAP ... # residues 32-45 only
rg_all:  GYRATION ATOMS=...
rg_nterm:GYRATION ATOMS=...
sasa:    SASA_HASEL  ATOMS=... TYPE=TOTAL
dom_sep: DISTANCE ATOMS=<N-domain COM>,<intermediate-domain COM>
PRINT ARG=q_all,q_nterm,q_loop,rg_all,rg_nterm,sasa,dom_sep FILE=COLVAR STRIDE=500
```

Build the contact maps from the folding units already computed:
`results/structure/folding_units.json` gives 6 spectral units and the per-residue contact order;
unit 2 (segments 1–4, 21–65, 80–89, 242–245) carries 17 of the 45 divergent positions including
all 14 loop residues, so `q_unit2` is worth adding as its own CV.

### (b) Umbrella sampling along domain separation and native-contact fraction

Gives the PMF of partial unfolding — the quantity that determines intermediate populations.

```bash
# pull the N-terminal GTPase domain away from the rest, 0.0 -> 2.5 nm
gmx grompp -f mdp/pull.mdp -c npt.gro -r npt.gro -p ${ID}.top -o pull.tpr
gmx mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg
python scripts_local/extract_windows.py --pullx pullx.xvg --spacing 0.1 --out windows/
for w in windows/*; do
  gmx grompp -f mdp/umbrella.mdp -c $w/start.gro -r $w/start.gro -p ${ID}.top -o $w/umb.tpr
  gmx mdrun -deffnm $w/umb -nb gpu -pme gpu
done
gmx wham -it tpr-files.dat -if pullf-files.dat -o pmf.xvg -hist hist.xvg -unit kCal -bsSteps 200
```

26 windows × 40 ns = 1.04 µs per system, ≈ 8 days on 1 A100. Check histogram overlap before
trusting the PMF; if adjacent windows do not overlap, add windows rather than extending.

Better if you have PLUMED comfort: 2-D well-tempered metadynamics on `(q_nterm, dom_sep)` with
`BIASFACTOR=15`, `SIGMA=0.02,0.05`, `HEIGHT=1.2`, reweighted afterwards. Same cost, and it gives
the 2-D surface rather than a 1-D cut.

### (c) High-temperature unfolding and quench — seeds only

```bash
# 10 x 20 ns at 450 K, then quench to 310 K for 20 ns each
for i in $(seq 1 10); do
  gmx grompp -f mdp/unfold450.mdp -c npt.gro -p ${ID}.top -o hot_$i.tpr
  gmx mdrun -deffnm hot_$i -nb gpu -pme gpu
done
```

These are **not** physical pathways — use them only as diverse starting points for Phase 3.

---

## Phase 3 — Markov state model: populations and lifetimes of the intermediates

This is what converts the ensemble into the numbers the kinetic model needs. The 50 000-draw
sensitivity sweep in `results/kinetics/` identified `k_degF` (clearance of the chaperone-free
quasi-native state) as the strongest discriminator of the phenotype, with a required shift of
**+0.925 log₁₀ units**. An MSM gives the state populations and mean first-passage times that this
can be checked against.

```bash
# seed 500 x 50 ns from clustered REST2 / quench structures (25 us aggregate)
python scripts_local/seed_trajectories.py --from rest2_${ID} --n 500 --ns 50 --out seeds/
# ... run them, then:
python scripts_local/build_msm.py \
  --trajs 'seeds/*/traj.xtc' --top ${ID}.pdb \
  --features backbone-torsions,q_unit2,dom_sep,loop_sasa \
  --tica-lag 10 --tica-dim 8 --n-clusters 400 --msm-lag 50 \
  --out msm_${ID}
```

Report, per isoform: number of metastable states (PCCA+), population of each, the native-contact
fraction and loop SASA of each, MFPT native → each intermediate and back, and the implied
lifetime of the most populated non-native state. **Validate with a Chapman–Kolmogorov test and an
implied-timescale plot before quoting anything** — an MSM that fails CK is worse than no MSM.

The comparison that matters: does TUBA8 populate a non-native, solvent-exposed intermediate more
than TUBA1A, and does the loop swap move that population?

---

## Phase 4 — does TBCB bind the intermediates differently?

The AF3 complexes model a **fully folded** substrate (α-tubulin chain pTM 0.91) and are suspect
for a further reason: the predicted TBCB interface covers **12 of the 18** α-tubulin residues
that β-tubulin uses in the mature heterodimer, burying 58% of the intradimer interface, with
**zero** overlap with the experimental TBCE interface in 9M1M. Either that is the real capping
mechanism or it is AF3 docking onto the largest hydrophobic patch. Phase 4 settles it.

### 4.1 Ensemble docking onto the MSM representatives

```bash
conda activate dock
# one HADDOCK run per representative intermediate (typically 5-10 states)
for state in msm_${ID}/rep_*.pdb; do
  haddock3 cfg/tbcb_alpha.cfg --extend-run runs/$(basename $state .pdb)
done
```

`cfg/tbcb_alpha.cfg`: `molecules = [tbcb.pdb, $state]`, restraints from the **conserved** part of
the AF3 interface only (the 49 residues shared between the TUBA8 and TUBA1A predictions — see
`results/af3/interface_aligned.json`), as ambiguous interaction restraints, with
`sampling_factor` raised so each state gets ≥ 1000 rigid-body models. Do **not** restrain to the
divergent positions: that would assume the answer.

Also run a restraint-free control (`CM` mode or AutoDock-Vina protein–protein) so you can tell
whether the restraints are driving the result.

### 4.2 MD refinement and MM/PBSA ranking

```bash
conda activate tubulin
# 5 x 20 ns per surviving pose, then
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ci index.ndx -cg 1 2 -ct traj.xtc -o FINAL.dat -eo per_frame.csv
```

`mmpbsa.in`: `startframe/endframe` over the last 15 ns, `interval=10`, `&gb igb=8` **and**
`&pb istrng=0.15, inp=2` so GB and PB are both reported, `&decomp idecomp=2, dec_verbose=0` for
per-residue decomposition. Quote the **spread across the 5 replicas**, not the per-frame standard
error, which is optimistic by an order of magnitude. MM/PBSA ranks poses; it does not give
affinities.

### 4.3 Interaction-energy and force decomposition

`scripts/interface_decomposition.py` in this repository does the exact pairwise MM decomposition
(Coulomb + Lennard-Jones per residue pair, per-residue totals, net force per residue, residual
force after minimisation). Run it on each refined complex:

```bash
python scripts/interface_decomposition.py --cif refined_complex.pdb --tag <name> --cutoff 12
```

---

## Phase 5 — rigorous free energies, in priority order

### 5.1 ΔΔG_fold(I114S) — highest return, run this first

A single-atom-scale change is the regime where alchemy converges. The thermodynamic cycle needs
both legs: the folded protein and an unfolded reference (the 5-residue fragment 112–116, already
built and validated in `gpu/systems/`).

```bash
python gpu/10_fep_mutation.py --system systems/a1a_solv.xml --mutation I114S \
  --legs folded --replicas 3 --lambdas 20 --ns-per-window 10 --out out/A1
python gpu/10_fep_mutation.py --system systems/fragment_I114.xml --mutation I114S \
  --legs reference --replicas 3 --lambdas 20 --ns-per-window 20 --out out/A2
python gpu/90_analyse_fep.py --job out/A1
```

Electrostatics and sterics decoupled across separate λ stages (a buried hydrophobic-to-polar
change has a large vdW component that single-stage schedules converge badly). Hamiltonian
replica exchange across the ladder. **Gates**: adjacent-window overlap ≥ 0.03 everywhere,
forward/reverse hysteresis within the replica spread, uncertainty from the spread across three
replicas. `gpu/90_analyse_fep.py` refuses to print a ΔΔG if the overlap gate fails.

≈ 40 GPU-h per leg. Also run I110V, G232S, A294S and A270V as controls (≈ 170 GPU-h total).

**Interpretation benchmark**: the kinetics require a `k_degF` shift of +0.925 dex. A folding
destabilisation of ≈ 1.3 kcal/mol corresponds to roughly a 9-fold change in the folded/unfolded
ratio at 310 K, so a ΔΔG_fold of 1–2 kcal/mol would make the mechanism quantitatively consistent;
anything below ~0.5 kcal/mol would rule I114S out.

### 5.2 Relative binding free energy, TUBA1A → TUBA8 at the interface

Only the tail positions differ at the predicted interface (440, 444, 446, 451), and the tail is
disordered in the middle (PAE 10–13 Å) while the terminal aromatic is anchored (PAE 3.6–3.8 Å).
So run the RBFE on the **terminal residue only**, in the complex and free:

```bash
python gpu/10_fep_mutation.py --system systems/tbcb_a1a_solv.xml --mutation Y451F \
  --legs complex free --replicas 3 --lambdas 12 --ns-per-window 10 --out out/C1
```

≈ 60 GPU-h. Expect a small number: in the 2HL3 crystal the tyrosine hydroxyl's nearest acceptor
is **3.68 Å**, outside hydrogen-bond range, and CAP-Gly specificity is canonically `-EEY/F`.

### 5.3 Absolute binding free energy, TBCB + α-tubulin — expensive, do last

Double decoupling with Boresch orientational restraints, or separated topologies. With a ~5200 Å²
interface this is a 2–4 µs commitment per system and the dominant error is the restraint
correction. Only attempt it if Phase 4 shows a pose difference worth quantifying; otherwise the
PMF route below is cheaper.

Alternative: PMF along centre-of-mass separation with orientational restraints, 40 windows ×
50 ns = 2 µs, plus the standard-state correction
`ΔG° = −ΔG_PMF − k_BT ln(V_site/V°)` — get the site-volume term right or the whole number is
meaningless.

### 5.4 The loop swap — not alchemically tractable

A 14-residue segment substitution needs the population comparison from Phases 2–3, not a single
ΔΔG. `gpu/20_rest2_loop.py` is set up for the REST2 side of it.

---

## Phase 6 — close the loop with the kinetic model

Every free energy and MSM rate above maps onto a parameter in
`scripts/maturation_model.py`. Re-run the sensitivity sweep with the measured values substituted:

```bash
python scripts/global_sensitivity_fast.py      # edit BASE with the measured values first
python scripts/selectivity_analysis.py
python scripts/regionality_sensitivity.py
```

The model's exact steady-state solver (`scripts/exact_steady.py`, validated to 10⁻¹³ against
stiff integration, 35–2000× faster) makes a 50 000-draw sweep take minutes, so this is free once
the physics is in hand.

The specific prediction to test: the phenotype region requires **kon_B +0.73 dex and
koff_B −0.25 dex**, i.e. **TUBA8 binding TBCB more tightly than TUBA1A** — the opposite of the
starting hypothesis, and currently unsupported by the predicted structures, whose interfaces are
superimposable with none of TUBA8's 38 core substitutions in contact. If Phase 4 finds no
difference either, the mechanism is not TBCB affinity and Phase 5.1 (`k_degF` via stability)
carries the explanation alone.

---

## Minimum viable programme

If GPU time is limited, this order maximises information per hour:

| order | job | cost | decides |
|---|---|---|---|
| 1 | ΔΔG_fold(I114S) + fragment leg (5.1) | ≈ 80 GPU-h | whether the one unbuffered substitution is the mechanism |
| 2 | REST2 on `a8` and `a1a` (2a) | ≈ 8 GPU-days | whether the intermediate ensembles differ at all |
| 3 | MSM from those (3) | CPU-heavy, ≈ 25 µs aggregate | populations and lifetimes, i.e. `k_degF` directly |
| 4 | ensemble docking + MM/PBSA (4.1, 4.2) | ≈ 3 GPU-days | whether TBCB prefers TUBA8 intermediates |
| 5 | Y451F RBFE (5.2) | ≈ 60 GPU-h | closes the EEY/EEF question quantitatively |
| 6 | loop-swap REST2 (2a on the swapped pair) | ≈ 8 GPU-days | the sharpest sequence/phenotype correlation |

Jobs 1 and 2 together answer the question the whole project has narrowed to. Everything before
them in this repository is sequence, structure and kinetics; none of it can substitute.

## What must not be reported from any of this without the stated gate

- an MM/PBSA number as an affinity or a K_d
- a ΔΔG whose λ-window overlap fell below 0.03, or whose forward and reverse runs disagree by
  more than the replica spread
- an MSM population or MFPT from a model that failed the Chapman–Kolmogorov test
- a docking score as an energy
- a PMF-derived ΔG° without the explicit site-volume standard-state correction
