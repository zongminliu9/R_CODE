# What the AF3 models can and cannot support, and where the divergence sits in the folding process

Prompted by two objections: that AF3 predicts end states which may simply be wrong, and that the
physiologically relevant species is a folding intermediate rather than native α-tubulin. Both are
correct, and acting on them changes the standing of several conclusions.

---

## 1. A concrete red flag in the AF3 complexes

AF3 was given full-length TBCB and full-length α-tubulin and returned a complex in which
α-tubulin is **fully folded** (chain pTM 0.91). TBCB does not act on native α-tubulin; it acts on
a quasi-native intermediate released from CCT. So the model answers "what would native
α-tubulin look like bound to TBCB", which may not be a physiological species at all. AF3's
confidence is calibrated on PDB-like end states, and a chaperone bound to a non-native substrate
is out of that distribution, so ipTM 0.83 does not license a mechanistic reading.

An independent check makes this concrete. Comparing the AF3 TBCB interface with the experimental
partner surfaces in 9M1M:

| pair | shared residues | Jaccard |
|---|---|---|
| TBCB (AF3) vs **TBCE** (9M1M) | **0** | 0.000 |
| TBCB (AF3) vs TBCD (9M1M) | 2 | 0.038 |
| **TBCB (AF3) vs β-tubulin** | **12** | **0.226** |

**The AF3 TBCB interface covers 12 of the 18 α-tubulin residues that β-tubulin uses in the mature
heterodimer, burying 393 Å² — 58% of the intradimer interface.** Two readings are consistent with
that and the data here cannot separate them:

- *Mechanistic*: TBCB caps the β-binding surface precisely to prevent premature heterodimer
  formation, so TBCB release must precede assembly. This is a sensible chaperone mechanism and
  would give the kinetic model's E→H step a physical basis.
- *Artefact*: the intradimer interface is the largest exposed hydrophobic patch on free
  α-tubulin, and AF3's training set contains thousands of αβ-tubulin structures. Docking a
  partner onto a genuine protein-binding surface belonging to a *different* partner is a known
  AF3 failure mode.

Also notable: TBCB's own chain pTM is **0.51**. The protein is substantially disordered on its
own, which is what one expects of a chaperone that must adapt to a non-native substrate — and it
is a direct warning that a rigid-body docked model is the wrong representation.

### The EM map cannot adjudicate it

EMD-2447 (α-tubulin + TBCE + TBCB) is the only experimental structure of a TBCB-containing
complex. It is 21 Å nominal on a 32³ grid with 4.66 Å voxels: **1 638 voxels for 136.8 kDa, about
1.3 voxels per residue.** Residue-level interface assignment needs better than 5 Å; this is four
times coarser.

What it does say: a single connected envelope of 165 757 Å³ with principal extents 85 × 83 × 50 Å
(Rg 37.2 Å), against the AF3 TBCB+α sub-complex at 83 × 48 × 41 Å (Rg 30.1 Å) — smaller in every
dimension, so it fits. That is a very weak statement.

Two honest caveats on this check. The volume bookkeeping ("space left over matches TBCE") is
**circular**: the contour threshold was chosen to enclose the total expected mass, so the
remainder necessarily matches the third component. And the reported 115 Å extent for TBCB alone is
an artefact of the principal-axis metric applied to a splayed, partly disordered chain, not a real
dimension. **The map provides essentially no constraint on the AF3 model.**

## 2. The intermediate-state objection, turned into a number

If the real substrate is partially folded, residues buried in the native state would be exposed
and could contact TBCB — which would reopen the kon_B hypothesis that §1 of report 03 had closed.
The test is how far each divergent position is from TBCB in the AF3 complexes (min heavy-atom
distance, mean over 5 models, sd ≤ 2.4 Å throughout):

| distance to TBCB | divergent positions |
|---|---|
| in contact (<4.5 Å) | 4 — all in the tail (440, 444, 446, 451) |
| 4.5–10 Å | 5 |
| 10–20 Å | 12 |
| **> 20 Å** | **22 (51.2%)** |

Core (non-tail) median distance is **21.1 Å**. The positions within 10 Å are only
I75V, T80A, I110V, M425L, V437T, V440F.

**The entire 32–45 loop sits 29–37 Å from TBCB** — the furthest cluster in the protein, on the
opposite face. Bringing it into contact needs a 30 Å displacement, which is domain-scale
unfolding, not the modest rearrangement a quasi-native intermediate implies. The
intermediate-state argument does not rescue the loop.

One convergence does stand out: **I110V at 9.0 Å and I114S at 13.1 Å** are the closest core
divergences to TBCB, and I114S is independently the single unbuffered destabilising substitution
(report 03, §4).

## 3. Where the divergence sits in the folding process

Rather than guessing at intermediate structures, three proxies that need no trajectory
(`scripts/folding_units.py`), on the 437 resolved residues of 9M1M chain a with 2 361 heavy-atom
contacts below 5 Å:

- **relative contact order** — mean sequence separation of a residue's contacts. Low means its
  environment is sequence-local, so it forms early and autonomously; high means distant parts of
  the chain must be in place first.
- **folding units** by spectral clustering of the contact graph (K = 6), giving compact
  self-contacting blocks, with the inter-unit boundary fraction per residue.
- **coevolutionary degree** from the DCA already computed, as the evolutionary signature of
  residues the rest of the fold depends on.

| metric | divergent set | rest | p |
|---|---|---|---|
| relative contact order | 0.0433 | 0.0675 | **p(low) = 0.0036** |
| unit-boundary fraction | 0.1463 | 0.1634 | 0.69 / 0.31 (n.s.) |
| coevolutionary degree | 14.38 | 12.99 | p(high) = 0.0224 |

| metric | **32–45 loop** | rest | p |
|---|---|---|---|
| relative contact order | **0.0185** | 0.0669 | **p(low) = 0.0000** |
| unit-boundary fraction | **0.0000** | 0.1672 | **p(low) = 0.0000** |
| coevolutionary degree | **17.75** | 12.96 | **p(high) = 0.0001** |

The loop has **the lowest contact order in the protein and zero contacts to any other folding
unit**. It is a sequence-local, self-contained element that forms early and independently, and it
does not participate in the late cooperative steps where a chaperone would be required. It is a
coevolutionary hub, but its couplings go to other TUBA8-divergent positions rather than to the
functional surfaces (report 03, §2).

Unit 2 (segments 1–4, 21–65, 80–89, 242–245) carries 17 of the divergent positions (24.3% of its
70 residues) including all 14 loop residues; no other unit exceeds 10.5%.

### A third buffering layer, and one convergence

Divergent positions sit at **low** contact order (p = 0.0036) — in the parts of the fold that form
early and locally. With the two constraints already found, TUBA8's substitution pattern is
buffered three ways:

| constraint | evidence |
|---|---|
| avoids buried positions | randomising positions *raises* the burial cost, z = −3.35 |
| preserves net charge | 45 divergent positions contribute **exactly 0** net charge change |
| avoids late-forming cooperative positions | relative contact order p(low) = 0.0036 |

**TUBA8's 45 substitutions are systematically placed where they do least damage, to both the
folding process and the final state.** The divergence is engineered to be tolerable.

Finally, a convergence worth recording. The chain-average relative contact order is **0.0653**,
low for a 437-residue globular protein (typically 0.08–0.15) — topologically this is an easy fold.
Yet tubulin is obligately chaperone-dependent. So tubulin's chaperone requirement is **not** a
topological folding-difficulty problem, which is exactly what the kinetic model found
independently: `k_fold` has no effect on isoform selectivity anywhere in a 50 000-draw prior.

## 4. Standing of the hypotheses after this round

| candidate | status | decisive evidence |
|---|---|---|
| TBCB interface difference | **no basis in the native-state model, and the model is suspect** | interfaces superimposable; but AF3 covers 58% of the β site and TBCB pTM is 0.51 |
| intermediate-state exposure rescuing kon_B | **excluded for the loop, open for 110–114** | loop is 29–37 Å away; I110V 9.0 Å, I114S 13.1 Å |
| long-range allostery | not supported | DCA 0.0 percentile; ANM agrees |
| C-terminal tail recognition | real contact, not differential | PAE 3.6 vs 3.8 Å; TUBA4A preserved with no aromatic |
| 32–45 loop acting on folding | **weakened at every test** | lowest contact order, zero inter-unit contacts, forms early and autonomously |
| **cumulative load + I114S** | **the strongest survivor** | total burial cost 3.4× the next isoform; I114S 2.4× the next single insult and 13 Å from TBCB |

The honest summary is that the biophysical hypotheses have been eliminated one by one, and the
substitution pattern turns out to be actively buffered rather than damaging. That makes the
remaining possibilities narrow and specific:

1. **Cumulative marginal destabilisation** despite per-position buffering — 17.2 cost units
   against 5.1 for TUBA4A. Marginal stability is additive, so many buffered insults can still sum
   to a short-lived chaperone-free state. This is the `k_degF` mechanism and is what the GPU
   package tests.
2. **I114S specifically** — the one unbuffered insult, and conveniently a single-atom-scale
   change where alchemical free energy converges.
3. **Something outside the folded-state/folding-pathway domain entirely.** The loop carries K40,
   and TUBA8 is the only α-tubulin that cannot be K40-acetylated. That is a post-translational
   rather than a folding difference, and nothing computed here would detect it.
