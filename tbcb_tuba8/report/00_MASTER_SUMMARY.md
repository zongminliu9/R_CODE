# TBCB-dependent TUBA8 stability — master summary

Final state of every claim, with the numbers that support it. The four working reports
(`01`–`04`) are kept in chronological order including conclusions later overturned; this document
supersedes them where they disagree.

**Before using any number: read `package/NOT_COMPUTED.md`.** There are no binding free energies,
no K_d values, no alchemical ΔΔG and no docking scores in this project. All 512 scalar results
carry `is_free_energy = False`.

---

## 1. The observation, and the one control that decides its interpretation

TBCB knockout markedly reduces TUBA8 protein in cerebellum; TUBA1A, TUBA1B and TUBA4A are
comparatively preserved, and TUBA8 itself is not reduced in hippocampus, prefrontal cortex or
striatum.

Proteomics: TUBA8 (Q9JJZ2) is validly quantified in **all 60 tissue TMT channels** across all
five regions. So TUBA8 is present where the knockout has no effect, and the regional restriction
is a mechanistic fact, not an expression floor.

That matters because the kinetic model's null test is decisive: with both isoforms given
**identical rate constants** and differing only in synthesis flux — scanned from 1× down to
0.01× — selectivity is **exactly 0.000000 at every ratio**. Mass-action competition for a shared
limiting chaperone pool produces no isoform selectivity whatsoever. **A TUBA8-specific rate
difference is required.**

## 2. TUBA8 is the extreme outlier of the family

Exact Needleman–Wunsch global alignment (BLOSUM62, affine −11/−1) over reviewed mouse sequences:

| pair | % identity | substitutions |
|---|---|---|
| TUBA8 vs TUBA1A | 90.42 | 43 |
| TUBA8 vs TUBA1B | 90.42 | 43 |
| TUBA8 vs TUBA1C | 90.20 | 44 |
| TUBA8 vs TUBA3A | 91.09 | 40 |
| TUBA8 vs TUBA4A | 90.18 | 44 |
| TUBA1A vs TUBA1B | 99.56 | 2 |
| TUBA3A vs TUBA4A (most divergent non-TUBA8 pair) | 94.64 | 24 |

TUBA8 sits at ~90% from *every* other isoform while the rest of the family spans 94.6–99.6%.
Against TUBA1A: **45 divergent canonical positions, 34 of them TUBA8-unique.** All four residues
named at the outset are real (K40→A, I110→V, A270→V, C-terminal Y→F), but they are 4 of 45, so
every analysis was run over the full set.

Divergence is concentrated in two places:

- **the H1–S2/K40 loop**: 12 substitutions in residues 32–45, i.e. 27% of the divergence in 3% of
  the sequence. Permutation null for the maximum 15-residue window density over 2×10⁴ draws:
  mean 4.79, 99th percentile 7, **maximum ever observed 9** — the observed 12 was never reached,
  **p = 5×10⁻⁵**. `PDGQMPSDKTIGGG` → `ADGTFGTQASKIND`, losing both prolines and two glycines.
  This segment is **strictly conserved in all five other mouse α-tubulins and replaced only in
  TUBA8** — a perfect correlation with the phenotype.
- **the C-terminal tail**: 6 divergent positions at 440–451, unresolved in every experimental
  structure (α-tubulin ends at 437 in 9M1M, 439 in 1JFF).

## 3. The substitution pattern is buffered in three independent ways

| constraint | measurement | p |
|---|---|---|
| avoids buried positions | keeping TUBA8's substitution identities but randomising their positions **raises** the burial cost: observed 17.23, null 22.27 ± 1.50 | z = −3.35, p(higher) = 0.9991 |
| preserves net charge | the 45 divergent positions contribute **exactly 0** net charge change (whole protein TUBA1A −24, TUBA8 −23) | — |
| avoids late-forming cooperative positions | relative contact order 0.0433 vs 0.0675 for the rest | p(low) = 0.0036 |
| enriched on exposed surface instead | 17 observed vs 11.16 expected | p = 0.0266 |

Per-substitution burial cost is unremarkable (TUBA8 0.442, against TUBA1C 0.436 and TUBA1B 0.577);
the large total simply reflects having 39 resolved substitutions instead of 2–13.

**TUBA8's 45 substitutions are placed where they do least damage to both the folding process and
the final state.** The divergence is engineered to be tolerable — which makes the phenotype
harder, not easier, to explain, and is the central tension of this project.

One substitution escapes the buffering:

| rank | substitution | burial | ΔFP (kcal/mol) | ΔVolume (Å³) | cost |
|---|---|---|---|---|---|
| **1** | **I114S** | **0.911** | **1.84** | **77.7** | **2.830** |
| 2 | I42K | 0.323 | 2.79 | −1.9 | 1.163 |
| 3 | G44N | 0.671 | 0.60 | −54.0 | 1.157 |
| 4 | G232S | 0.985 | 0.04 | −28.9 | 0.998 |
| 5 | A294S | 0.855 | 0.35 | −0.4 | 0.985 |

**I114S is 2.4× the next-largest single insult**: a 91%-buried isoleucine replaced by serine,
losing 1.84 kcal/mol of hydrophobic burial and 78 Å³ of volume — it both opens a cavity and
buries a hydroxyl. It sits in a local cluster (I110V, I114S, A126T, Q128A, T130S in 21 residues)
and is 13.1 Å from TBCB, the second-closest core divergence.

Composition note: **12 of 43 substitutions go to Ser or Thr (27.9%)** against a background S+T
frequency of 11.3% — expected 4.9, **binomial p = 0.0023**; to Ser alone 7 vs 2.0, p = 0.0035;
net +6 hydroxyl side chains. They are **not** spatially clustered (mean pairwise Cα 29.54 Å vs
null 27.96 ± 2.44, p = 0.74) and **not** a distinguishable coevolutionary sub-module
(z = +0.97, p = 0.17).

## 4. The substitutions are not independent drift

Mean-field DCA on the 44 191-sequence α-tubulin MSA, validated first: **60.7% contact precision
on the top L/4 long-range pairs (|i−j| ≥ 24), 62.9× enrichment over a 0.97% baseline**.
Meff/L = 2.42.

Within-set coupling of the 45 divergent positions, against three nulls:

| null | mean | sd | z | p |
|---|---|---|---|---|
| uniform random positions | −0.0158 | 0.0098 | +18.94 | 5×10⁻⁵ |
| **entropy-matched positions** | 0.0697 | 0.0173 | **+5.76** | **5×10⁻⁵** |
| coupling-matrix label shuffle | 0.0055 | 0.0110 | +14.98 | 2×10⁻⁴ |

The entropy bias is real and large (divergent-set per-site entropy 1.425 vs 0.694) and absorbs
most of the naive z. The effect survives it, and survives deleting both spatial clusters so that
only the **26 scattered single-site substitutions** remain: **z = +4.44, p = 4×10⁻⁴**.

Evolution treats TUBA8's substitutions as a coupled module. What they are coupled *to* is each
other — including long-range pairs 42↔80 (29 Å), 45↔159 (31 Å) and **41↔437 (54 Å)**, the K40
loop coevolving with the C-terminal region — and not to the functional surfaces.

## 5. Every structural mechanism for differential TBCB handling has failed

### 5.1 Folded-core cofactor interfaces — excluded

Permutation over the 437 resolved α-tubulin residues of 9M1M (10⁵ draws):

| structural class | observed | expected | p |
|---|---|---|---|
| cofactor interface (ΔSASA > 5 Å²) | **0** | 4.91 | p_deplete = **0.0041** |
| within 8 Å of TBCE | 3 | 8.30 | p_deplete = **0.0187** |
| nucleotide pocket (<5 Å) | **0** | 2.59 | 0.0617 |
| α–β intradimer interface | **0** | 1.61 | 0.1783 |

Scope limit: 9M1M contains **no TBCB**.

### 5.2 The TBCB interface itself — no basis for a kon_B or koff_B difference

From the AF3 complexes (ipTM 0.830 TUBA8 / 0.8160 TUBA1A, every TUBA8 model ≥ every TUBA1A
model; buried SASA 5250.6 ± 90.7 vs 5206.5 ± 158.0 Å², a difference **within the scatter**):

- TBCB wraps α-tubulin with all three of its parts — ubiquitin-like domain 18–19 contacts,
  linker 21, **CAP-Gly 18–21 including Lys203 and Asn204, the GKNDG carboxylate readers**.
- Compared in canonical numbering, the two interfaces are **superimposable**: 49 of 52/56
  residues shared; divergent positions at the interface **identical in both**,
  [440, 444, 446, 451]; permutation p_enrich 0.79, p_deplete 0.39 — at chance.
- **4 of 7 tail divergent positions contact TBCB; 0 of 38 core ones do.**

### 5.3 The intermediate-state rescue — excluded for the loop, open for 110–114

Distance from each divergent position to TBCB (min heavy-atom, mean over 5 models, sd ≤ 2.4 Å):

| distance | divergent positions |
|---|---|
| in contact (<4.5 Å) | 4, all in the tail |
| 4.5–10 Å | 5 |
| 10–20 Å | 12 |
| **> 20 Å** | **22 (51.2%)** |

Core median **21.1 Å**. **The entire 32–45 loop is 29–37 Å away**, on the opposite face: reaching
TBCB needs domain-scale unfolding, not the modest rearrangement a quasi-native intermediate
implies. Closest core divergences: **I110V 9.0 Å, I114S 13.1 Å**.

### 5.4 Long-range allostery — not supported, by the test that could have found it

The elastic-network result alone could not carry this: an ANM is harmonic and geometric, blind to
entropy-driven and side-chain-mediated coupling, and with a symmetric covariance its PRS matrix
is symmetric so effectiveness and sensitivity are identical and carry no directional information.
It gave only non-significance (p = 0.16–0.22, 78th–84th percentile of matched windows).

The coevolutionary test is emphatic. Restricted to pairs **far apart in space** (Cβ–Cβ > 12 Å,
because coupling between residues that touch is a contact, not allostery), null = every other
contiguous 14-residue window:

| target surface | loop coupling | null mean | **percentile** |
|---|---|---|---|
| **TBCB interface (AF3)** | −0.1607 | −0.0223 | **0.0** |
| α–β interface | −0.1681 | −0.0186 | 0.7 |
| nucleotide pocket | −0.1508 | −0.0255 | 4.4 |
| late-cofactor interface | −0.0649 | −0.0147 | 32.0 |

The loop is **less coupled to the TBCB interface than any other 14-residue window in the
protein**. None of the 15 strongest effectors onto the cofactor patch is a divergent position.

### 5.5 C-terminal tail recognition — a real contact, but not differential

The tail **is** the TBCB contact (this overturns report 01). PAE shows a tethered anchor:

| | TUBA8 | TUBA1A |
|---|---|---|
| tail middle | 9.1–11.0 Å (440–443) | 10.6–12.9 Å (441–445) |
| **terminal aromatic → CAP-Gly** | **F449: 3.60 ± 0.06 Å** | **Y451: 3.84 ± 0.21 Å** |

Three independent lines say EEF does not weaken engagement:

1. **2HL3 crystal geometry.** Recognition is carried by the free carboxylate — Phe88 backbone
   amide at **2.64 Å**, the GKNDG Asn69 ND2 at **2.74 Å** — plus ring packing on Phe52 (3.89 Å).
   The tyrosine hydroxyl's nearest acceptor is **3.68 Å**, outside hydrogen-bond range. It is not
   a specificity determinant, consistent with CAP-Gly specificity being canonically `-EEY/F`.
2. **AF3 confidence.** TUBA8's Phe terminus is placed marginally *more* confidently than
   TUBA1A's Tyr.
3. **The TUBA4A control.** TUBA4A ends `…EDEDEGEE` — a glutamate, **no terminal aromatic at
   all**, unreadable by any CAP-Gly domain — yet it is preserved. With a ~5200 Å² interface built
   almost entirely from the conserved core, the terminal anchor is one contact among many. This
   argument uses sequence and crystal geometry only, no force field.

Common-geometry energetics agree in sign but point the *wrong way* for the phenotype:
Y→F costs **+2.91 kcal/mol at TBCB** (losing a 2.80 Å hydroxyl–Phe187 O contact, itself probably
a GB over-contact given the crystal's 3.68 Å, so this is an upper bound) and **+0.00 at TBCE**,
where the hydroxyl has no polar neighbour within 4.5 Å. That would make TUBA8 *less* TBCB-dependent.

### 5.6 The 32–45 loop acting on folding — weakened at every test

| metric | 32–45 loop | rest | p |
|---|---|---|---|
| relative contact order | **0.0185** | 0.0669 | p(low) = **0.0000** |
| contacts to other folding units | **0.0000** | 0.1672 | p(low) = **0.0000** |
| coevolutionary degree | 17.75 | 12.96 | p(high) = 0.0001 |

The loop has the **lowest contact order in the protein and zero inter-unit contacts**: a
sequence-local, self-contained element that forms early and autonomously, and therefore does not
participate in the late cooperative steps where a chaperone is required. Two cheap proxies for
loop disorder were tried and both are invalidated by their own controls (pLDDT control spread
10.72 against a TUBA8 deviation of 0.76; relative SASA 0.379–0.526 across three
identical-sequence isoforms).

## 6. What the kinetics require, and the tension that remains

50 000-draw sweep (49 728 usable), broad prior: every baseline rate ±1 decade, every
TUBA8/TUBA1A ratio ±1.5 decades, demand / TBCB_tot / TBCE_tot over 3 decades. Phenotype
(sel < −1.5, FC₈ < 0.5, FC₁ₐ > 0.75) reached in **4.39%** of draws.

Controlling parameters (Spearman ρ):

| parameter | ρ | p |
|---|---|---|
| **kon_Edir** — TBCB-independent TBCE capture | **+0.445** | 0 |
| **k_degF** — clearance of unprotected quasi-native | **−0.293** | 0 |
| **kon_B** — TBCB capture | **−0.255** | 0 |
| koff_B | +0.105 | 5×10⁻¹²³ |
| **k_fold** | **−0.007** | **0.11, n.s.** |

Phenotype-region median shifts (log₁₀): **k_degF +0.925, kon_Edir −1.079, kon_B +0.731, koff_B −0.247, TBCB_tot +0.387, demand −0.280**. TBCE_tot does not discriminate (p = 0.059).

Three readings matter:

- **The phenotype needs a conjunction**, not one parameter. This corrects the one-at-a-time
  reading: `k_degF` saturates near 25% loss alone but is the strongest discriminator in
  combination, because a `kon_Edir` deficit keeps TUBA8 in the unprotected state longer, which is
  what makes F's clearance rate decisive.
- **Folding rate is irrelevant everywhere in the prior.** Independently, the chain-average
  relative contact order is **0.0653** — low for 437 residues, i.e. topologically an easy fold —
  yet tubulin is obligately chaperone-dependent. Tubulin's chaperone requirement is not a
  folding-difficulty problem. Two methods, same conclusion.
- **The regional variable is the TBCB-to-substrate ratio, not TBCE capacity**, and **lower**
  synthesis flux favours the phenotype (ρ = −0.059; region median −0.280), the opposite of a
  naive "cerebellum has the highest tubulin demand" expectation. The |sel| ridge tracks
  demand ≈ TBCE_tot.
- **11.9% of phenotype draws need no kon_Edir deficit at all** (median k_degF +1.038,
  kon_B +1.048, koff_B −0.494), so the observation does not uniquely imply a missing bypass route.

**The tension.** The kinetics require kon_B ↑ and koff_B ↓ for TUBA8, i.e. TUBA8 binding TBCB
*better*. The predicted interface offers no basis for that: the interfaces are superimposable and
none of the 38 core substitutions touches TBCB. This is reported unresolved. The AF3 caveat in §7
is the most likely way out, but it is not demonstrated.

## 7. How far the AF3 models can be trusted

The models supply the interface that no experiment maps, and they must be read with three
reservations:

1. **They model a fully folded substrate** (α-tubulin chain pTM 0.91), whereas TBCB acts on a
   quasi-native intermediate from CCT. The complex answers "what would *native* α-tubulin look
   like bound to TBCB", which may not be a physiological species.
2. **A concrete red flag.** The predicted TBCB interface covers **12 of the 18** α-tubulin
   residues that β-tubulin uses in the mature heterodimer, burying **393 Å² = 58%** of the
   intradimer interface, with **zero** overlap with the experimental TBCE interface in 9M1M
   (Jaccard 0.000). Either TBCB caps the β site to block premature dimerisation — a sensible
   chaperone mechanism that would give the E→H step a physical basis — or AF3 docked onto the
   largest exposed hydrophobic patch, which is a known failure mode. **The data here cannot
   separate these.**
3. **TBCB's own chain pTM is 0.51**: it is substantially disordered alone, so a rigid-body docked
   model is the wrong representation regardless.

EMD-2447 cannot adjudicate: 21 Å nominal, 32³ grid, 4.66 Å voxels, **~1.3 voxels per residue**
against the <5 Å needed for interface assignment. The envelope is a single lobe of 85 × 83 × 50 Å
(Rg 37.2 Å) and the AF3 sub-complex (83 × 48 × 41 Å, Rg 30.1 Å) fits inside it — a very weak
statement. Two of our own checks on it are flawed and flagged as such: the volume bookkeeping is
**circular** (the contour was chosen to match total mass) and the 115 Å extent reported for TBCB
alone is a principal-axis artefact on a splayed chain.

## 8. Standing of every hypothesis

| hypothesis | verdict | decisive evidence |
|---|---|---|
| low TUBA8 expression explains the regional pattern | **excluded** | proteomics: all 60 TMT channels; null model gives selectivity exactly 0.000000 |
| folded-core cofactor interface difference | **excluded** | 0 divergent residues at the 9M1M cofactor interface, p = 0.0041 |
| TBCB binding affinity difference (kon_B / koff_B) | **no structural basis found** | interfaces superimposable; 0 of 38 core substitutions contact TBCB |
| intermediate-state exposure rescuing the above | **excluded for the loop, open for 110–114** | loop 29–37 Å from TBCB; I110V 9.0 Å, I114S 13.1 Å |
| long-range allostery from the loop | **not supported** | DCA 0.0 percentile; ANM agrees |
| C-terminal EEY→EEF recognition | **excluded as a differential** | TUBA4A preserved with no aromatic; hydroxyl at 3.68 Å makes no H-bond; energetics point the wrong way |
| higher folding barrier | **excluded** | k_fold ρ = −0.007 (p = 0.11) across 49 728 draws; contact order says the fold is easy |
| slower TBCB→TBCE handoff | **excluded, wrong sign** | reducing k_tr *decreases* TBCB dependence |
| 32–45 loop acting on the folding pathway | **weakened at every test** | lowest contact order, zero inter-unit contacts, forms early and autonomously |
| **cumulative marginal destabilisation (k_degF)** | **strongest survivor** | total burial load 3.4× the next isoform; the only mechanism the sweep supports in conjunction |
| **I114S specifically** | **sharpest single-residue lead** | 2.4× the next insult, 91% buried, 13 Å from TBCB |
| loss of K40 acetylation (post-translational, not folding) | **untested here** | TUBA8 is the only α-tubulin that cannot be K40-acetylated; nothing computed would detect this |

## 9. What would settle it

**Two predictions testable with no further simulation:**

1. **TUBA8 should bind TBCB more tightly than TUBA1A**, not less (kon_B +0.73 dex, koff_B −0.25 dex in the phenotype region). This is the opposite of the starting hypothesis and
   is directly measurable.
2. **The 32–45 loop should be swappable.** Replacing TUBA8's 32–45 with the conserved
   `PDGQMPSDKTIGGG` should rescue TUBA8 in a TBCB-null background, and the reciprocal swap into
   TUBA1A should confer dependence. The perfect sequence/phenotype correlation makes this the
   sharpest available test — though §5.6 now predicts it will fail, which is itself informative.

**Two free-energy calculations, inputs built and validated in `gpu/`, not run** (no GPU reachable
from the analysis environment; measured throughput in `package/benchmarks.csv`):

- **ΔΔG_fold(I114S)** by alchemical FEP, 20 λ windows with electrostatics and sterics decoupled,
  three replicas, MBAR with an overlap floor of 0.03 and forward/reverse hysteresis reported.
  ≈170 GPU-h including the I110V, G232S, A294S and A270V controls. A destabilisation of
  ≳1.5 kcal/mol would correspond to roughly a 12-fold shift in the folded/unfolded ratio at
  310 K, which is the order the kinetics require (+0.925 dex in k_degF).
- **REST2 comparison of 32–45 loop local unfolding**, TUBA8 against the sequence-swapped control.
  ≈400 GPU-h.

TBCB binding free energy is deliberately **not** prioritised: with superimposable interfaces and
no core substitution in contact, there is no structural hypothesis for a ΔΔG_bind to test.

**And one line outside biophysics.** Every structural and energetic analysis says the
substitutions are engineered to be harmless. If that is right, TUBA8's special status may not lie
in its own physics but in a regulatory site it has **lost**: the 32–45 loop carries K40, and
TUBA8 is the only α-tubulin that cannot be K40-acetylated. That is a post-translational
difference, it would explain a regional phenotype if cerebellar α-tubulin acetylation flux
differs, and **nothing computed in this project would detect it.**

> **Correction (applied after a numerical cross-check).** An earlier version of this
> file quoted the phenotype-region *medians* (k_degF 0.879, kon_Edir −1.014, kon_B 0.689,
> TBCB_tot 0.379, k_tr 0.236, koff_B −0.227) while labelling them *shifts*. The shift is
> the phenotype-region median minus the median over the remaining draws, and the correct
> values are +0.925, −1.079, +0.731, +0.387, +0.241 and −0.247 log10 units. The machine-
> readable tables in `package/` always carried the correct shifts; only the prose was
> mislabelled. Conclusions are unchanged — the ranking and signs are identical.
