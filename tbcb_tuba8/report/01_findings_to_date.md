# TBCB-dependent TUBA8 stability — computational findings to date

Workspace: `/home/user/tbcb-tuba8`. Every number below is produced by a script in `scripts/`
and written to `results/`; nothing here is quoted from memory or from the literature.

---

## 0. Executive summary of what the computations actually say

Three premises of the starting plan survive, three do not.

**Survives.** TUBA8 really is the extreme outlier of the mouse α-tubulin family, and all four
pre-named residue differences (A40, I110V, A270V, C-terminal F rather than Y) are real.

**Does not survive.** (a) The divergent residues are *statistically depleted* at every
functional interface in the late cofactor complex, including the TBCE contact surface, and
*enriched* on the free solvent-exposed surface. So a model in which TUBA8's TBCB dependence
comes from altered chaperone-contact chemistry in the folded core is not supported by the
structural distribution of its substitutions. (b) In the kinetic model, a higher folding
barrier (`k_fold`) produces **exactly zero** selectivity, so "TUBA8 folds more slowly" cannot
by itself be the mechanism. (c) A slower TBCB→TBCE handoff (`k_tr`) changes selectivity in the
**wrong direction**.

**What the computations point to instead.** TUBA8's divergence is concentrated in two
intrinsically flexible, solvent-exposed places — a near-wholesale replacement of the H1–S2/K40
loop (residues 32–45) and the disordered C-terminal tail — and the only kinetic parameter that
can generate the observed magnitude of selective loss is the rate of **TBCB-independent capture
by TBCE**. Together these say the question to ask is not "does TBCB bind TUBA8 more tightly"
but "can TUBA8 reach TBCE without TBCB".

---

## 1. Isoform divergence (`scripts/align_isoforms.py`)

Exact Needleman–Wunsch global alignment (BLOSUM62, affine gaps, open −11 / extend −1) over
reviewed mouse UniProt sequences. α-tubulins are >90% identical, so the global optimum is
unambiguous and no heuristic aligner is needed.

| pair | % identity | substitutions |
|---|---|---|
| TUBA8 vs TUBA1A | 90.42 | 43 |
| TUBA8 vs TUBA1B | 90.42 | 43 |
| TUBA8 vs TUBA1C | 90.20 | 44 |
| TUBA8 vs TUBA3A | 91.09 | 40 |
| TUBA8 vs TUBA4A | 90.18 | 44 |
| TUBA1A vs TUBA1B | 99.56 | 2 |
| TUBA3A vs TUBA4A (most divergent non-TUBA8 pair) | 94.64 | 24 |

TUBA8 sits ~90% from *every* other isoform while the rest of the family spans 94.6–99.6%.
Against TUBA1A there are **45 divergent canonical positions, 34 of which are TUBA8-unique**
(shared with no other mouse α-tubulin). The starting plan named four; the unbiased count is 45,
so all downstream analysis was run over the full set rather than the pre-selected four.

Verification of the pre-named positions (canonical TUBA1A/1JFF numbering): K40→A **confirmed**,
I110→V **confirmed**, A270→V **confirmed**, C-terminus `...GEEEGEEY` → `...EENEGEEF`
**confirmed** (TUBA4A ends `...EDEDEGEE`, lacking the aromatic altogether).

## 2. Where the divergence sits (`scripts/divergence_vs_interfaces.py`, `scripts/enrichment_test.py`)

Structural context computed from **9M1M** (TBCD–TBCE–TBCC–ARL2–αβ-tubulin, cryo-EM 2.21 Å).
Chains were assigned by sequence identity to the UniProt references rather than hard-coded, and
resolve as C=TBCC, D=TBCD, E=TBCE, G=ARL2, a=α-tubulin, b=β-tubulin. The α-site holds GTP·Mg²⁺;
the β-site holds GDP·AlF₃.

Per-residue SASA (Shrake–Rupley, 200 points, Tien et al. 2013 Gly-X-Gly normalisation) was
computed in three contexts — α alone, α+β, and α+β+all cofactors — to separate intradimer
burial from cofactor burial. Permutation null: 10⁵ random position sets of the same size drawn
from the 437 resolved α-tubulin residues.

| structural class | class size | observed | expected | ratio | p(enrich) | p(deplete) |
|---|---|---|---|---|---|---|
| cofactor interface (ΔSASA>5 Å²) | 55 | **0** | 4.91 | 0.00 | 1.000 | **0.0041** |
| within 8 Å of TBCE | 93 | 3 | 8.30 | 0.36 | 0.995 | **0.0187** |
| nucleotide pocket (<5 Å) | 29 | **0** | 2.59 | 0.00 | 1.000 | 0.0617 |
| α–β intradimer interface | 18 | **0** | 1.61 | 0.00 | 1.000 | 0.1783 |
| buried core (relSASA<0.10) | 178 | 11 | 15.88 | 0.69 | 0.969 | 0.0653 |
| exposed surface (relSASA>0.35) | 125 | **17** | 11.16 | 1.52 | **0.0266** | 0.989 |

Not one divergent residue touches a cofactor in this complex, and the TBCE neighbourhood is
significantly depleted. The substitutions avoid the nucleotide pocket and the intradimer
interface entirely and accumulate on the free surface.

**Important scope limit.** 9M1M is a *late* maturation complex and **contains no TBCB**. These
numbers therefore constrain the TBCD/TBCE/TBCC/ARL2 surfaces of the late complex, not the
TBCB–α-tubulin interface, for which no experimental structure exists. The depletion result
weakens the "altered cofactor-contact chemistry" model but cannot by itself exclude a TBCB-
specific effect.

### The one dramatic exception: residues 32–45

Sliding 15-residue window over the full canonical 1–451:

| window | substitutions | positions |
|---|---|---|
| **32–46** | **12** | 32, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45 |
| 31–45 | 12 | same |
| 36–50 | 11 | 36–45, 50 |

Null distribution of the maximum window density for 45 randomly placed substitutions
(2×10⁴ draws): mean 4.79, 99th percentile 7, **maximum ever observed 9**. The observed value of
12 was never reached, giving **p = 5×10⁻⁵**.

So 12 of 45 differences (27%) fall in 3% of the sequence. `DKTIGGG` → `QASKIND` across 39–45 is
not a set of point mutations, it is effectively a **replacement of the H1–S2 / K40 acetylation
loop**. This is the single most distinctive feature of TUBA8 by a wide margin.

### The C-terminal tail is disordered in every available structure

Six divergent positions (440, 441, 442, 444, 446, 451) lie beyond the last residue resolved in
either reference structure — α-tubulin ends at 437 in 9M1M and at 439 in 1JFF. The EEY→EEF
difference therefore sits in an intrinsically disordered tail and **cannot** be analysed as a
static contact; it has to be treated as a conformational-ensemble plus binding problem.

Both loci of TUBA8 distinctiveness — the 32–45 loop and the C-terminal tail — are flexible and
solvent-exposed. Neither is a folded-core chaperone contact.

## 3. Multistate maturation kinetics (`scripts/maturation_model.py`, `scripts/selectivity_analysis.py`)

States per isoform: U (post-CCT, partly folded), F (quasi-native, cofactor-free), B (TBCB-bound),
E (TBCE-bound), H (mature soluble heterodimer), M (assembled/retained), D (cleared, absorbing).
TBCB and TBCE are **shared, conserved pools**, so isoforms compete:
`TBCB_tot = TBCB_free + Σ B_i`, `TBCE_tot = TBCE_free + Σ E_i`. Observable = F+B+E+H+M.
Selectivity `sel = log₂FC(TUBA8) − log₂FC(TUBA1A)`, `FC = KO/WT`. Steady states verified
stationary to a relative drift of ~10⁻¹⁸.

### The null model is rejected

With *identical* rate constants for both isoforms and TUBA8 differing only in synthesis flux —
scanned from 1× down to 0.01× of TUBA1A — selectivity is **0.0000 at every ratio**. Mass-action
competition for a shared limiting chaperone pool produces no isoform selectivity whatsoever;
both isoforms fall by the same few percent. **Low TUBA8 expression cannot explain the phenotype.**
A genuine TUBA8-specific rate difference is required.

### One mechanism at a time

| mechanism varied in TUBA8 only | strongest selectivity reachable | verdict |
|---|---|---|
| **kon_Edir** — TBCB-independent TBCE capture | **sel = −3.30** (FC₈ = 0.10, FC₁ₐ = 0.988) | **sufficient; does not saturate** |
| kon_B — TBCB capture rate | −0.46 (FC₈ = 0.73) | too weak, saturates |
| k_degF — clearance of unprotected F | −0.32 (FC₈ = 0.75) | too weak, saturates |
| k_fold — completion of folding | **0.000** | **no effect at all** |
| k_tr — TBCB→TBCE handoff | **+0.03** | **wrong sign** |
| k_asm — heterodimer completion | 0.000 | no effect |

Only a deficient TBCB-independent route to TBCE can generate selectivity of the observed
magnitude. "Marginal stability → faster clearance" (k_degF) and "higher folding barrier"
(k_fold) are quantitatively inadequate as standalone explanations — the first saturates near
25% loss, the second does nothing.

### Regional specificity must be a flux variable, and the model predicts a ridge, not a gradient

TUBA8's sequence is identical in cerebellum and cortex, so the regional restriction of the
phenotype has to come from a flux or abundance variable. Giving TUBA8 a single difference
(kon_Edir reduced 30×) and scanning synthesis demand:

| demand | FC_TUBA8 | FC_TUBA1A | sel |
|---|---|---|---|
| 0.03 | 0.259 | 0.910 | −1.81 |
| 0.3 | 0.215 | 0.886 | −2.04 |
| **1.0** | 0.222 | 0.985 | **−2.15** |
| 3 | 0.633 | 1.006 | −0.67 |
| 10 | 0.949 | 1.001 | −0.08 |
| 100 | 0.999 | 1.000 | −0.001 |

Selectivity **vanishes at high demand**, the opposite of a naive "cerebellum has the highest
tubulin demand" argument: once substrate flux saturates TBCE, the bottleneck moves downstream to
chaperone turnover and the delivery-route difference washes out.

The 2-D scan over demand × TBCE abundance shows the maximum of |sel| tracking the diagonal
**demand ≈ TBCE_tot** (peak at demand 1 for TBCE 1, at 3 for TBCE 3, at 10 for TBCE 10). The
prediction is therefore that the phenotype appears where α-tubulin synthesis flux is *matched*
to TBCE capacity; regions with large chaperone excess **and** regions with flux far above
capacity are both protected.

---

## 4. Computational feasibility of the originally proposed simulation programme

Measured on this container (4 × Xeon 2.10 GHz, 15 GB RAM, **no GPU**; OpenMM 8.6.1, CPU platform):

| system | atoms | conditions | throughput |
|---|---|---|---|
| CAP-Gly + EEY peptide (2HL3) | 1 227 | GBn2 implicit, 2.0 nm cutoff, 4 fs | 6.2 ns/day |
| CAP-Gly + EEY peptide | 1 227 | GBn2, 1.0 nm cutoff, 4 fs | 20.3 ns/day |
| CAP-Gly + EEY peptide | 1 227 | GBn2, 1.2 nm, HMR, 5 fs | 17.1 ns/day |
| mouse TBCB CAP-Gly (1WHG) | 1 650 | GBn2, 2.0 nm, 4 fs | 3.6 ns/day |
| α-tubulin monomer (1JFF) | 6 375 | GBn2, 1.0 nm, HMR, 5 fs | **3.9 ns/day** |

Single-thread throughput for the 1 227-atom system is 5.15 ns/day, so four λ windows run in
parallel give ~20 ns/day aggregate.

Consequences, with the *implicit-solvent shortcut already applied* (the plan specified explicit
ff19SB/OPC, which is a further 10–30× slower):

- Simulation A1 as written — 4 isoforms × 3 replicas × 20 ns of α-tubulin monomer — needs
  **62 days** of this machine.
- Protein–protein PMF (Simulation B3): 30–50 windows × 20–100 ns ≈ 1–5 μs → decades.
- Per-mutation FEP (Simulation C2): ~400 ns of λ sampling → months per mutation.
- Simulations D, E and G involve systems up to ~2 500 residues (9M1M), i.e. larger still.

The programme as specified is on the order of 10³–10⁴ GPU-days. **The binding free energies,
PMFs and FEP ΔΔG values it calls for are not obtainable here, and no number of that kind should
be reported from this container.** They are quoted above only as resource estimates.

### What this hardware can do rigorously

1. Comparative sequence/structure determinant analysis — done, sections 1–2.
2. Interface burial and permutation statistics — done, section 2.
3. Elastic-network / normal-mode analysis and perturbation-response scanning: exact linear
   algebra on a ~1 300×1 300 matrix, seconds, no sampling error. The right substitute for the
   "conformational dynamics" aim.
4. Minimisation-based interaction energies and local strain — cheap, but they are energies, not
   free energies, and must be labelled as such.
5. Short (1–5 ns) MD of the CAP-Gly + tail peptide system for contact occupancy — hours.
6. The full multistate kinetic/thermodynamic model with sensitivity analysis — done, section 3.
7. A GPU-ready package (prepared inputs, force-field settings, λ schedules, protocols) so the
   MD/FEP programme can be executed unchanged on appropriate hardware.

---

## 5. Data provenance

- 49 reviewed UniProt entries (mouse + human), `data/seq/`, index in `data/meta/`.
- 30 experimental structures, `data/pdb/`; 8 AlphaFold v6 models, `data/afdb/`.
- Key structures located: **9M1M** 2.21 Å late cofactor complex; **9M1L/9M1J** related states;
  **2HL3** p150Glued CAP-Gly with a bound **EEY tripeptide** (the binding-mode template for the
  EEY/EEF question); **1WHG** mouse TBCB CAP-Gly and **1V6E** mouse TBCB ubiquitin-like domain
  (both domains of the mouse protein, experimentally); **4B6M** TBCB CAP-Gly at 1.59 Å;
  **4ICV** TBCE ubiquitin-like domain.
- No experimental structure of a TBCB–α-tubulin complex exists; that interface is unmapped.

## 6. The control the model says matters most

Because the null model is rejected only on the assumption that TUBA8 is *expressed* in the
unaffected regions, the decisive experimental control is whether TUBA8 protein is present at a
measurable level in hippocampus, cortex and striatum in wild type. If it is near zero there, the
absence of a knockout effect in those regions is an expression floor rather than a mechanism,
and no biophysical explanation is required for it.
