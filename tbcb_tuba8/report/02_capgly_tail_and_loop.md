# Testing the two candidate determinants: the C-terminal tail and the H1–S2/K40 loop

Follow-on to `01_findings_to_date.md` §2, which localised TUBA8's divergence to two
intrinsically flexible, solvent-exposed regions, and §8, which showed the kinetic phenotype
requires a conjunction of ↑kon_B, ↓koff_B, ↓↓kon_Edir and ↑↑k_degF.

**Outcome: the C-terminal tail hypothesis is refuted. The 32–45 loop survives as the only
candidate, by elimination, and with a perfect phenotype correlation — but the computation that
would actually test it needs sampling this hardware cannot provide.**

---

## 1. Both cofactors read the tail, and they differ at the reading residue

The kinetic model needs an asymmetry between TBCB capture (`kon_B`) and direct TBCE capture
(`kon_Edir`). Both cofactors carry CAP-Gly domains, so both read the α-tubulin
`-EEY/F-COO⁻` terminus — which makes differential reading a structurally concrete candidate.

Motif search over the reviewed sequences and the coordinates (`scripts/find_tbce_capgly.py`):

| protein | CAP-Gly motif | position | carboxylate-reading residue |
|---|---|---|---|
| p150Glued (2HL3, the reference complex) | **GKNDG** | 67–71 | **Asn69** |
| mouse / human TBCB | **GKNDG** | 202 | **Asn204** |
| mouse / human **TBCE** | **GK<u>H</u>DG** | 47 | **His49** |

TBCB keeps the canonical asparagine; TBCE has a histidine in its place. Mouse TBCB's CAP-Gly
spans Q9D1E6 residues **134–231** (established by mapping the 1WHG NESG construct, auth 8 =
residue 134, after removing the GSSGSSG / SGPSSG tags). TBCE's spans residues **4–121** and is
resolved in 9M1M chain E at 2.21 Å.

In 9M1M only **one** TBCE residue (108) lies within 12 Å of α-tubulin 430–437, so the CAP-Gly
domain sits where the tail emerges but the tail is not captured in that state.

## 2. The experimental complex says the hydroxyl is not a specificity determinant

Contacts of the bound EEY tripeptide in 2HL3 (`scripts/validate_transplant.py`); 0 clashes,
138 atom pairs < 4.5 Å, so this is the physical baseline:

| Tyr atom | nearest domain atom | distance |
|---|---|---|
| **OXT** (free carboxylate) | Phe88 **N** (backbone amide) | **2.64 Å** |
| **O** | Asn69 **ND2** (the GKNDG asparagine) | **2.74 Å** |
| **OH** (the group Phe lacks) | Phe52 **O** (backbone carbonyl) | **3.68 Å** |
| CZ (ring) | Phe52 CD2 | 3.89 Å |

Recognition is carried by the free carboxylate, read at 2.6–2.7 Å by the Phe88 backbone amide
and the GKNDG asparagine, plus ring packing against Phe52. The tyrosine hydroxyl's nearest
potential acceptor is **3.68 Å** away — outside normal hydrogen-bond range. It makes no good
hydrogen bond, which is consistent with the canonical description of CAP-Gly specificity as
`-EEY/F`, where phenylalanine is an accepted ligand rather than a rejected one.

## 3. Transplanted models, and an energetic estimate that had to be discarded

The EEY binding mode was transplanted into each cofactor groove by iterative Kabsch
superposition on sequence-aligned Cα with outlier rejection (`scripts/build_four_complexes.py`):

| model | template | identity to p150 | Kabsch RMSD | clashes < 2 Å | pairs < 4.5 Å |
|---|---|---|---|---|---|
| TBCB CAP-Gly | AlphaFold Q9D1E6 134–231 | 52.0% | **0.47 Å** / 56 | 4 | 176 |
| TBCE CAP-Gly | 9M1M chain E 4–121 | 34.2% | 1.10 Å / 50 | 3 | 126 |
| T. brucei TBCB (check) | 4B6M, 1.59 Å | 40.3% | 0.87 Å / 55 | 6 | 134 |
| mouse TBCB (check) | 1WHG NMR | 52.0% | 1.00 Å / 56 | 32 | 492 |

The contact sets map correctly onto the native pattern (W192↔W57, K203/N204↔K68/N69,
F223↔F88, V224↔V89 for TBCB; W36↔W57, K48/H49↔K68/N69, F69↔F88 for TBCE), so the transplants
are geometrically sound. 1WHG is overpacked because its NMR side chains clash, and 4B6M
reproduces the native pattern almost exactly (134 pairs against 138).

**A first energetic attempt was discarded as unreliable** (`scripts/interaction_energies.py`).
Minimising the EEY and EEF complexes independently gave ΔΔE_int(vac) = −27.5 kcal/mol for TBCB
from deleting a single hydroxyl, which is not physical: the two runs converged to different
local minima. Compounding this, ΔE_int(vac) ≈ −300 kcal/mol is nearly cancelled by
ΔG_solv ≈ +265 kcal/mol for a peptide of net charge −3, so the result is a small difference of
large numbers and a 1% error in either term moves it by ±3 kcal/mol. Single-structure MM/GBSA
cannot resolve a ΔΔG of this size. The numbers are retained in the results directory only as a
record of a rejected approach.

### The common-geometry probe

Replacing it with a protocol that removes the dominant error source
(`scripts/hydroxyl_probe.py`): minimise the EEY complex once, then build PHE in place by
replacing the hydroxyl oxygen with the HZ it becomes, placed along the CZ–OH vector at 1.08 Å,
leaving every other atom where it was.

| groove | polar atoms within 4.5 Å of the hydroxyl | Y→F ΔΔE(GB) |
|---|---|---|
| **TBCB** | **Phe187 O at 2.80 Å**, plus Lys188 aliphatic packing (2.70–2.75 Å) | **+2.91 kcal/mol** |
| **TBCE** | **none** — only Val32 methyls at 3.45–4.23 Å | **+0.00 kcal/mol** |

Caveat on the magnitude: the crystal has this contact at 3.68 Å and restrained minimisation in
GB implicit solvent pulled it to 2.80 Å. With no explicit water to compete, GB models
systematically over-form direct hydrogen bonds, so +2.91 kcal/mol should be read as an upper
bound on the hydroxyl's contribution at TBCB.

### This asymmetry points the wrong way

The phenotype requires TUBA8 to be captured **well** by TBCB and **poorly** by TBCE. The
calculation gives the opposite: the EEF terminus is disadvantaged by up to ~2.9 kcal/mol at
TBCB and exactly unchanged at TBCE. On this evidence TUBA8 should be *less* TBCB-dependent,
not more.

## 4. The decisive argument is in the sequences, and it needs no force field

| isoform | C-terminal 16 residues | last | CAP-Gly readable? | TBCB KO |
|---|---|---|---|---|
| TUBA1A | GVDSVEGEGEEEGEEY | Y | yes | preserved |
| TUBA1B | GVDSVEGEGEEEGEEY | Y | yes | preserved |
| TUBA1C | EVGADSAEGDDEGEEY | Y | yes | preserved |
| TUBA3A | VGVDSVEAEAEEGEEY | Y | yes | preserved |
| **TUBA4A** | EEVGIDSYEDEDEGE**E** | **E** | **no terminal aromatic at all** | **preserved** |
| **TUBA8** | EVGTDSFEEENEGEE**F** | **F** | yes (F is an accepted ligand) | **LOST** |

TUBA4A ends in a glutamate and cannot be read by a CAP-Gly domain at all, yet it is preserved.
TUBA8 ends in a phenylalanine, which CAP-Gly accepts, yet it is lost. If differential CAP-Gly
reading of the tail set TBCB dependence, TUBA4A would be the most affected isoform and TUBA8
among the least. The observed ordering is the reverse.

**The C-terminal tail is not the determinant.** This conclusion rests on sequence and on the
experimental 2HL3 geometry, not on any force field, minimisation protocol or implicit-solvent
model, so it is more robust than the energetics in §3 — which happen to agree with it.

## 5. What remains: the H1–S2/K40 loop, and a perfect correlation

| isoform | residues 32–45 | Gly | Pro | TBCB KO |
|---|---|---|---|---|
| TUBA1A | `PDGQMPSDKTIGGG` | 4 | 2 | preserved |
| TUBA1B | `PDGQMPSDKTIGGG` | 4 | 2 | preserved |
| TUBA1C | identical to TUBA1A | 4 | 2 | preserved |
| TUBA3A | identical to TUBA1A | 4 | 2 | preserved |
| TUBA4A | `PDGQMPSDKTIGGG` | 4 | 2 | preserved |
| **TUBA8** | **`ADGTFGTQASKIND`** | **2** | **0** | **LOST** |

The segment is **strictly conserved across all five other mouse α-tubulins and replaced only in
TUBA8** — a perfect correlation with the phenotype, and precisely the property the C-terminal
tail lacks. TUBA8 loses both prolines and two glycines and gains Phe36 plus a set of polar and
charged side chains, turning a Pro/Gly-rich, conformationally restricted loop into an ordinary
one bearing side chains.

### Two cheap tests of the loop, both without discriminating power

Honest reporting of two metrics that did **not** work, each invalidated by its own control.

**Hydrophobic exposure** (`scripts/loop32_45_properties.py`), relative loop SASA in the
AlphaFold models: TUBA1A 0.379, TUBA1B 0.526, TUBA4A 0.474, TUBA8 0.456. The three isoforms with
*identical* loop sequence span 0.379–0.526, a spread of 0.15 that exceeds any TUBA8 deviation.
The metric is dominated by model noise in a flexible loop.

**AlphaFold pLDDT** (`scripts/plddt_loop.py`):

| isoform | loop pLDDT | rest of chain | Δ |
|---|---|---|---|
| TUBA1A | 76.81 | 91.48 | −14.67 |
| TUBA1B | 66.09 | 92.23 | −26.14 |
| TUBA4A | 69.89 | 92.82 | −22.93 |
| **TUBA8 (mouse)** | **70.17** | 92.48 | −22.31 |
| TUBA8 (human) | 69.94 | 91.90 | −21.97 |

Control spread over the three identical-sequence isoforms is 66.09–76.81 (range 10.72), and
TUBA8 sits 0.76 below their mean — **well within the noise**. Mouse and human TUBA8 agree to
0.23, so this is not single-prediction scatter. The loop is a low-confidence region in *every*
α-tubulin (≈66–77 against ≈92 elsewhere): it is intrinsically flexible in all of them and TUBA8
is unexceptional. AlphaFold does not predict TUBA8's loop to be more disordered.

## 6. Where this leaves the investigation

Of the four candidate loci, three are excluded:

| locus | status | evidence |
|---|---|---|
| folded-core cofactor interfaces | **excluded** | divergence depleted, p = 0.004; zero contacts in 9M1M |
| long-range allosteric coupling | **excluded** | no significant coupling, p = 0.16–0.22; no divergent residue among the top effectors |
| C-terminal tail / CAP-Gly recognition | **excluded** | TUBA4A control; hydroxyl makes no H-bond at 3.68 Å; energetics point the wrong way |
| **H1–S2/K40 loop (32–45)** | **the only survivor** | 12 substitutions, p = 5×10⁻⁵; strictly conserved in all five other isoforms; perfect phenotype correlation |

The kinetic model says what the loop would have to do: not slow folding (`k_fold` has no effect
on selectivity anywhere in the prior) but **shorten the lifetime of the chaperone-free
quasi-native state** (`k_degF`, the strongest discriminator of the phenotype region, +0.879 dex).

Testing that means comparing the free-energy landscape and clearance-relevant exposure of the
chaperone-free quasi-native state between TUBA8 and TUBA1A — the free-energy problem of the
original Simulation A, which at 3.9 ns/day for an α-tubulin monomer this container cannot
deliver (§4 of the first report). The two cheap proxies attempted here are, by their own
controls, incapable of substituting for it.

So the useful outcome is a sharpened and much smaller question than the one we started with:
**does replacing the Pro/Gly-rich H1–S2 loop raise the clearance rate of chaperone-free
α-tubulin?** That is a single comparison on a defined 14-residue substitution, it has a clear
computational protocol, and it needs GPU time rather than more analysis.

### Two predictions that can be checked without any simulation

1. **TUBA8 should bind TBCB more tightly than TUBA1A does**, not less (§8 of the first report:
   the phenotype region requires kon_B +0.69 dex and koff_B −0.23 dex). This is the opposite of
   the starting hypothesis and is directly measurable.
2. **The 32–45 loop should be transplantable.** Replacing TUBA8's 32–45 with the conserved
   `PDGQMPSDKTIGGG` should rescue TUBA8 levels in a TBCB-null background, and the reciprocal
   swap into TUBA1A should confer TBCB dependence. The perfect sequence/phenotype correlation
   in §5 makes this the sharpest available test, and it requires no free-energy calculation.
