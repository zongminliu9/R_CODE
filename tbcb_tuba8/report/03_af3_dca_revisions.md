# Revisions forced by the AF3 complexes, the DCA, and the proteomics result

This supersedes parts of `02_capgly_tail_and_loop.md`. Three earlier conclusions change, one is
strengthened by a much better test, and one new single-residue lead appears.

---

## 0. What changed

| earlier claim | status now | why |
|---|---|---|
| "no experimental or predicted structure maps the TBCB interface" | **obsolete** | AF3 complexes supplied, ipTM 0.83/0.82, interface ≈ 5200 Å² |
| "the C-terminal tail is disordered, so it is not a TBCB contact" | **wrong** | the tail **is** the contact; the terminal aromatic is docked in the CAP-Gly groove at PAE 3.6–3.8 Å |
| "the tail hypothesis is excluded" | **overstated → reframed** | EEF does not weaken engagement, but the tail is genuinely part of the interface |
| "long-range allosteric coupling is excluded" | **overstated → now properly supported** | ANM alone could not exclude it; DCA puts the loop at the **0.0 percentile** of coupling to the TBCB interface |
| "the regional phenotype may be an expression floor" | **excluded** | TUBA8 quantified in all 60 TMT channels across all five regions |
| — | **new lead** | **I114S**, a buried Ile→Ser, is the only substitution in TUBA8 that is not buffered |

The proteomics answer matters: TUBA8 is present in hippocampus, cortex, striatum and spinal cord,
so the absence of a knockout effect there is a real mechanistic fact that needs explaining, not an
artefact of TUBA8 being absent. The §3 null-model result — that equal rate constants with
different synthesis flux give selectivity **exactly 0.0000** — therefore stands unweakened, and a
TUBA8-specific rate difference is required.

## 1. The TBCB–α-tubulin interface (AF3, 5 models per complex)

Confidence. ipTM is the interface score; pTM is inflated by the two folded chains and must not
substitute for it.

| complex | ipTM (5 models) | mean | chain pTM (TBCB, α) | clashes |
|---|---|---|---|---|
| TBCB + TUBA8 | 0.83 ×5 | **0.8300** | 0.51, 0.91 | 0 |
| TBCB + TUBA1A | 0.82, 0.82, 0.82, 0.81, 0.81 | 0.8160 | 0.50, 0.91 | 0 |

Every TUBA8 model scores ≥ every TUBA1A model. ipTM is a confidence measure, not an affinity, so
the 0.014 gap is not a ΔΔG — but its consistency across five models is worth recording because
its direction matches what the kinetic model demanded (TUBA8 binding TBCB at least as well).

**TBCB wraps α-tubulin with all three of its parts**, giving an unusually large interface
(≈5200 Å² buried, against 1200–2000 Å² for a typical complex):

| TBCB region | residues contacting α-tubulin |
|---|---|
| ubiquitin-like domain (1–99) | 18–19 |
| linker (100–133) | 21 |
| CAP-Gly domain (134–231) | 18–21, **including Lys203 and Asn204**, the GKNDG carboxylate readers |

### The interfaces are superimposable, and that is a problem for the model

Compared in canonical TUBA1A numbering (the first pass compared raw residue numbers and
manufactured false differences, since TUBA8 is 449 aa and TUBA1A 451):

- TUBA8 interface 52 residues, TUBA1A 56, **49 shared**
- divergent positions at the interface: **[440, 444, 446, 451] in both** — identical
- permutation test over canonical 1–451: observed 4, expected 5.19, **p_enrich 0.79, p_deplete 0.39** — at chance, neither enriched nor depleted (unlike the late 9M1M complex, where divergence was depleted at p = 0.004)

And the decisive split:

| divergent positions | at the TBCB interface |
|---|---|
| in the tail (437–451), 7 total | **4** |
| in the core (<437), 38 total | **0** |

Not one of TUBA8's 38 core substitutions touches TBCB. **The predicted interface provides
essentially no structural basis for a kon_B or koff_B difference between TUBA8 and TUBA1A**,
which is a genuine tension with the phenotype region's requirement of kon_B +0.69 dex and
koff_B −0.23 dex (first report, §8). That tension is reported, not resolved.

### The tail is a tethered anchor, and PAE shows it

PAE profile along the tail, minimum over TBCB residues, mean ± sd over 5 models:

| | TUBA8 | TUBA1A |
|---|---|---|
| middle of tail | 9.1–11.0 Å (440–443) | 10.6–12.9 Å (441–445) |
| **terminal aromatic → CAP-Gly** | **F449: 3.60 ± 0.06 Å** | **Y451: 3.84 ± 0.21 Å** |
| tail mean → CAP-Gly | 8.63 Å | 9.11 Å |

Pinned at both ends, genuinely disordered in between: the canonical picture of a flexible tether
with one specific anchor. The terminal `-EEY/F` is confidently docked (PAE < 4 Å) while the
intervening residues are not placed at all (PAE 10–13 Å).

TUBA8's Phe terminus is placed *marginally more* confidently than TUBA1A's Tyr. So the
EEY→EEF change does not weaken CAP-Gly engagement — now supported by two independent routes: the
2HL3 crystal geometry (the hydroxyl's nearest acceptor is 3.68 Å, outside H-bond range) and the
AF3 confidence profile.

### Resolving the TUBA4A paradox

TUBA4A has no terminal aromatic at all yet is preserved. With a ≈5200 Å² interface built almost
entirely from the conserved core, the terminal anchor is one contact among many, so losing it
costs comparatively little. That dissolves the paradox and simultaneously demotes the tail: the
terminal aromatic is probably not rate-limiting for TBCB capture by any isoform.

## 2. Allostery, tested properly

The earlier elastic-network result could not carry the weight I put on it. An ANM is harmonic and
purely geometric — blind to entropy-driven coupling, side-chain-mediated communication and any
coupling that is functional rather than mechanical. A null ANM result does not exclude long-range
coupling.

Direct coupling analysis asks a categorically different question: which residue pairs have
co-varied over evolution *more than transitive chains of other correlations explain*. Inverting
the covariance matrix is what removes those indirect chains.

**Setup and validation** (`scripts/dca2.py`): 44 191 match-state sequences from the AF3 paired
MSA, chunked reweighting at 80% identity (the full 44 191² identity matrix is 7.8 GB, so it is
never materialised), Meff = 1091.3, Meff/L = 2.42.

Contact precision against 9M1M chain a (Cβ–Cβ < 8 Å):

| separation | random baseline | top L/4 | enrichment |
|---|---|---|---|
| \|i−j\| ≥ 6 | 1.06% | 48.2% | 45.5× |
| \|i−j\| ≥ 12 | 0.98% | 54.5% | 55.5× |
| **\|i−j\| ≥ 24** | 0.97% | **60.7%** | **62.9×** |

60.7% precision on the top 112 long-range pairs is within the normal working range for
mean-field DCA, so the couplings carry real information despite the modest Meff.

### The 32–45 loop is maximally decoupled from the functional surfaces

Restricting to pairs that are **far apart in space** (Cβ–Cβ > 12 Å) and ≥6 apart in sequence,
because a coupling between residues that touch is a contact, not allostery. Null: every other
contiguous 14-residue window.

| target surface | loop coupling | null mean | **percentile** |
|---|---|---|---|
| **TBCB interface (AF3)** | −0.1607 | −0.0223 | **0.0** |
| α–β interface | −0.1681 | −0.0186 | 0.7 |
| nucleotide pocket | −0.1508 | −0.0255 | 4.4 |
| late-cofactor interface | −0.0649 | −0.0147 | 32.0 |

The loop is *less* coupled to the TBCB interface than **any** other 14-residue window in the
protein. Two methodologically independent probes — harmonic mechanics and evolutionary coupling —
now agree, and the second is emphatic rather than merely non-significant. A long-range allosteric
route from this loop to the cofactor machinery is not supported.

Its strongest genuinely long-range couplings go instead to *other TUBA8-divergent positions*:
42↔80 (29 Å), 45↔159 (31 Å), **41↔437 (54 Å)** — the K40 loop coevolving with the C-terminal
region, i.e. the two loci of TUBA8 divergence are coupled to each other.

### TUBA8's substitutions are a coevolving module — and this survives the control that mattered

Mean APC-corrected coupling within the 45-position divergent set, against three nulls:

| null model | mean | sd | z | p |
|---|---|---|---|---|
| uniform random positions | −0.0158 | 0.0098 | +18.94 | 5×10⁻⁵ |
| **entropy-matched positions** | 0.0697 | 0.0173 | **+5.76** | **5×10⁻⁵** |
| coupling-matrix label shuffle | 0.0055 | 0.0110 | +14.98 | 2×10⁻⁴ |

The entropy bias is real and large — divergent-set per-site entropy is 1.425 against 0.694
elsewhere, and matching it absorbs most of the naive z. **The effect survives**: after matching
conservation, and even after deleting both spatial clusters so that only the 26 scattered
single-site substitutions remain, within-set coupling is still enriched (**z = +4.44,
p = 4×10⁻⁴**).

So the scattered surface substitutions are **not independent drift**. Evolution treats them as a
coupled set. What they are coupled *to* is each other, not the functional surfaces.

## 3. The serine shift is real; its simplest interpretations are not

| | observed | expected | p |
|---|---|---|---|
| substitutions to Ser or Thr | **12 / 43 (27.9%)** | 4.9 (background S+T = 11.3%) | **0.0023** |
| substitutions to Ser alone | **7** | 2.0 (background S = 4.7%) | **0.0035** |

Net +6 hydroxyl side chains at divergent positions:
`Q35T S38T T41S N50T I114S A126T T130S G232S N293S A294S N300S V437T`.

Two follow-ups are negative, and both are informative:

- **Not spatially clustered**: mean pairwise Cα distance 29.54 Å against a null of 27.96 ± 2.44 Å (p = 0.74). They are scattered over the fold, so they are not building one patch or one composite site.
- **Not a distinguishable sub-module**: coupling among them is 0.2411 against random subsets *of the same divergent set* at 0.1627 ± 0.0808 (z = +0.97, p = 0.17). They belong to the general module rather than forming their own.

Two compositional constraints are worth recording:

- **Net charge is preserved almost exactly.** Whole protein: TUBA1A −24, TUBA8 −23. The 45 divergent positions contribute **exactly 0** net charge change. Forty-five substitutions with zero net charge change is a strong constraint signature for a protein whose lattice, MAP and motor interactions are electrostatic.
- Aromatic count is essentially unchanged (43 → 44).

## 4. Buried-substitution load: TUBA8 is buffered, with one exception

The kinetic model's surviving requirement (`k_degF`, faster clearance of the chaperone-free state)
is a folded-state *stability* property, not an interface property, so the relevant quantity is how
much buried packing each isoform's substitutions disturb. Proxy — and it is a proxy, not a free
energy — scored identically for every isoform against the same reference, with burial taken from
9M1M: `burial × hydrophobicity lost + 0.6 × burial × |ΔV|/100 + 0.8 × burial × newly-buried-polar`.

| isoform | resolved substitutions | total cost | **mean per substitution** |
|---|---|---|---|
| **TUBA8** | 39 | **17.23** | **0.442** |
| TUBA4A | 13 | 5.08 | 0.390 |
| TUBA3A | 8 | 2.32 | 0.290 |
| TUBA1C | 3 | 1.31 | 0.436 |
| TUBA1B | 2 | 1.15 | 0.577 |

TUBA8's total load is 3.4× the next isoform's, but its **per-substitution cost is unremarkable**
(0.442, between TUBA1C's 0.436 and TUBA1B's 0.577). The large total simply reflects having 39
substitutions instead of 2–13.

And the placement test is decisively negative in the opposite direction: keeping TUBA8's
substitution identities but randomising their positions gives 22.27 ± 1.50 against an observed
17.23 — **z = −3.35, p(higher) = 0.9991**. TUBA8's substitutions sit at positions *less buried
than random*. Evolution has buffered them. This is the thermodynamic restatement of the earlier
finding that the divergence is enriched on the exposed surface (p = 0.027).

Both things are true and should be stated together: the per-position insult is buffered, while the
cumulative load is still the largest in the family by 3.4×, and marginal stability is additive.

### The one unbuffered substitution: I114S

| rank | substitution | burial | ΔFP (kcal/mol) | ΔVolume (Å³) | cost |
|---|---|---|---|---|---|
| **1** | **I114S** | **0.911** | **1.84** | **77.7** | **2.830** |
| 2 | I42K | 0.323 | 2.79 | −1.9 | 1.163 |
| 3 | G44N | 0.671 | 0.60 | −54.0 | 1.157 |
| 4 | G232S | 0.985 | 0.04 | −28.9 | 0.998 |
| 5 | A294S | 0.855 | 0.35 | −0.4 | 0.985 |

I114S is 2.4× the next-largest single insult: a buried isoleucine replaced by serine, losing
1.84 kcal/mol of hydrophobic burial and 78 Å³ of side-chain volume at a position that is 91%
buried. It both opens a cavity and buries a hydroxyl. Substitutions of that kind typically cost
1.5–3 kcal/mol of folding stability.

It also sits in a local cluster: residues 110–130 carry five substitutions (I110V, I114S, A126T,
Q128A, T130S) in 21 residues, and two of the five — I110V and the segment around it — were among
the residues named at the outset.

**I114S is now the sharpest single-residue hypothesis in the project**, and unlike the loop swap
it is a single-atom-scale change, which is exactly the regime where an alchemical free-energy
calculation converges well.

## 5. Where the investigation stands

| candidate | status | best evidence |
|---|---|---|
| folded-core cofactor interfaces | excluded | 0/38 core divergent positions contact TBCB in the AF3 complexes; depleted in 9M1M (p = 0.004) |
| long-range allostery from the loop | **not supported** (strong) | DCA: 0.0 percentile coupling to the TBCB interface; ANM agrees |
| C-terminal tail recognition | **is a real contact, but not differential** | PAE 3.6 vs 3.8 Å; hydroxyl makes no H-bond; TUBA4A preserved without any aromatic |
| TBCB binding affinity difference | **no structural basis found** | interfaces superimposable, 49/52 shared, divergent-at-interface identical |
| **folded-state stability (k_degF)** | **the surviving mechanism** | largest cumulative buried load (3.4×); **I114S** unbuffered at 2.83 |
| 32–45 loop | **still correlates perfectly**, mechanism unclear | strictly conserved in all five other isoforms, replaced only in TUBA8 |

The two surviving candidates are now **I114S** and the **32–45 loop**, and they are testable by
the same quantity: the folding free energy and the clearance propensity of the chaperone-free
α-tubulin monomer. Both require free-energy sampling of a 450-residue protein, which this
container cannot do at 3.9 ns/day. The GPU package in `gpu/` is set up for exactly these two
calculations.
