# Data provenance

Every input, with accession and source. Structures and AlphaFold models are gitignored because
they are large and byte-for-byte regenerable by `scripts/fetch_sequences.py` and
`scripts/fetch_structures.py`.

## Sequences — UniProt, reviewed entries only

Retrieved by `scripts/fetch_sequences.py` via the UniProt REST API with
`(gene_exact:X) AND (organism_id:…) AND (reviewed:true)`.

| key | accession | gene | organism | length | protein |
|---|---|---|---|---|---|
| `human_ARL2` | P36404 | ARL2 | human | 184 | ADP-ribosylation factor-like protein 2 |
| `human_CCT2` | P78371 | CCT2 | human | 535 | T-complex protein 1 subunit beta |
| `human_CCT5` | P48643 | CCT5 | human | 541 | T-complex protein 1 subunit epsilon |
| `human_PFDN1` | O60925 | PFDN1 | human | 122 | Prefoldin subunit 1 |
| `human_PFDN2` | Q9UHV9 | PFDN2 | human | 154 | Prefoldin subunit 2 |
| `human_TBCA` | O75347 | TBCA | human | 108 | Tubulin-specific chaperone A |
| `human_TBCB` | Q99426 | TBCB | human | 244 | Tubulin-folding cofactor B |
| `human_TBCC` | Q15814 | TBCC | human | 346 | Tubulin-specific chaperone C |
| `human_TBCD` | Q9BTW9 | TBCD | human | 1192 | Tubulin-specific chaperone D |
| `human_TBCE` | Q15813 | TBCE | human | 527 | Tubulin-specific chaperone E |
| `human_TCP1` | P17987 | TCP1 | human | 556 | T-complex protein 1 subunit alpha |
| `human_TTL` | Q8NG68 | TTL | human | 377 | Tubulin--tyrosine ligase |
| `human_TUBA1A` | Q71U36 | TUBA1A | human | 451 | Tubulin alpha-1A chain |
| `human_TUBA1B` | P68363 | TUBA1B | human | 451 | Tubulin alpha-1B chain |
| `human_TUBA1C` | Q9BQE3 | TUBA1C | human | 449 | Tubulin alpha-1C chain |
| `human_TUBA4A` | P68366 | TUBA4A | human | 448 | Tubulin alpha-4A chain |
| `human_TUBA8` | Q9NY65 | TUBA8 | human | 449 | Tubulin alpha-8 chain |
| `human_TUBB` | P07437 | TUBB | human | 444 | Tubulin beta chain |
| `human_TUBB2A` | Q13885 | TUBB2A | human | 445 | Tubulin beta-2A chain |
| `human_TUBB2B` | Q9BVA1 | TUBB2B | human | 445 | Tubulin beta-2B chain |
| `human_TUBB3` | Q13509 | TUBB3 | human | 450 | Tubulin beta-3 chain |
| `human_TUBB4A` | P04350 | TUBB4A | human | 444 | Tubulin beta-4A chain |
| `human_TUBB4B` | P68371 | TUBB4B | human | 445 | Tubulin beta-4B chain |
| `human_VASH1` | Q7L8A9 | VASH1 | human | 365 | Tubulinyl-Tyr carboxypeptidase 1 |
| `mouse_Arl2` | Q9D0J4 | Arl2 | mouse | 184 | ADP-ribosylation factor-like protein 2 |
| `mouse_Cct2` | P80314 | Cct2 | mouse | 535 | T-complex protein 1 subunit beta |
| `mouse_Cct5` | P80316 | Cct5 | mouse | 541 | T-complex protein 1 subunit epsilon |
| `mouse_Pfdn1` | Q9CWM4 | Pfdn1 | mouse | 122 | Prefoldin subunit 1 |
| `mouse_Pfdn2` | O70591 | Pfdn2 | mouse | 154 | Prefoldin subunit 2 |
| `mouse_Tbca` | P48428 | Tbca | mouse | 108 | Tubulin-specific chaperone A |
| `mouse_Tbcb` | Q9D1E6 | Tbcb | mouse | 244 | Tubulin-folding cofactor B |
| `mouse_Tbcc` | Q8VCN9 | Tbcc | mouse | 341 | Tubulin-specific chaperone C |
| `mouse_Tbcd` | Q8BYA0 | Tbcd | mouse | 1196 | Tubulin-specific chaperone D |
| `mouse_Tbce` | Q8CIV8 | Tbce | mouse | 524 | Tubulin-specific chaperone E |
| `mouse_Tcp1` | P11983 | Tcp1 | mouse | 556 | T-complex protein 1 subunit alpha |
| `mouse_Ttl` | P38585 | Ttl | mouse | 377 | Tubulin--tyrosine ligase |
| `mouse_Tuba1a` | P68369 | Tuba1a | mouse | 451 | Tubulin alpha-1A chain |
| `mouse_Tuba1b` | P05213 | Tuba1b | mouse | 451 | Tubulin alpha-1B chain |
| `mouse_Tuba1c` | P68373 | Tuba1c | mouse | 449 | Tubulin alpha-1C chain |
| `mouse_Tuba3a` | P05214 | Tuba3a | mouse | 450 | Tubulin alpha-3 chain |
| `mouse_Tuba4a` | P68368 | Tuba4a | mouse | 448 | Tubulin alpha-4A chain |
| `mouse_Tuba8` | Q9JJZ2 | Tuba8 | mouse | 449 | Tubulin alpha-8 chain |
| `mouse_Tubb2a` | Q7TMM9 | Tubb2a | mouse | 445 | Tubulin beta-2A chain |
| `mouse_Tubb2b` | Q9CWF2 | Tubb2b | mouse | 445 | Tubulin beta-2B chain |
| `mouse_Tubb3` | Q9ERD7 | Tubb3 | mouse | 450 | Tubulin beta-3 chain |
| `mouse_Tubb4a` | Q9D6F9 | Tubb4a | mouse | 444 | Tubulin beta-4A chain |
| `mouse_Tubb4b` | P68372 | Tubb4b | mouse | 445 | Tubulin beta-4B chain |
| `mouse_Tubb5` | P99024 | Tubb5 | mouse | 444 | Tubulin beta-5 chain |
| `mouse_Vash1` | Q8C1W1 | Vash1 | mouse | 375 | Tubulinyl-Tyr carboxypeptidase 1 |

Total: **49 entries**. Primary reference for TUBA8 is mouse Q9JJZ2 as specified.

## Experimental structures — RCSB PDB

| PDB | description | method | resolution | notes |
|---|---|---|---|---|
| 9M1M | TBC-DEC-Arl2-ab-tubulin, GDP | EM | 2.21 A | chain assignment by sequence: C=TBCC, D=TBCD, E=TBCE, G=ARL2, a=alpha (TUBA1B, 1-437), b=beta (TUBB5, 1-431); alpha site GTP+Mg, beta site GDP+AlF3 |
| 9M1L | TBC-DE-Arl2-ab-tubulin, GTP | EM | 2.55 A |  |
| 9M1J | TBC-DE-Arl2-b-tubulin | EM | 2.14 A |  |
| 9EEB | TBCD/E/G + tubulin dimer | EM | 3.7 A |  |
| 9EDR | TBCD/E/G/C + tubulin, TBCC N-term | EM | 3.8 A |  |
| 1WHG | mouse TBCB CAP-Gly | NMR | - | NESG construct; GSSGSSG and SGPSSG tags removed; core auth 8-105 maps to Q9D1E6 residues 134-231 |
| 1V6E | mouse TBCB N-terminal ubiquitin-like domain | NMR | - |  |
| 4B6M | T. brucei TBCB CAP-Gly | X-ray | 1.59 A | motif GKGDG at 195-199 |
| 4B6W | T. brucei TBCB architecture | X-ray | 2.35 A |  |
| 1T0Y | TBCB ubiquitin-like domain | NMR | - |  |
| 4ICV | human TBCE ubiquitin-like domain | X-ray | 1.45 A |  |
| 2HL3 | p150Glued CAP-Gly (A49M) + EEY tripeptide | X-ray | - | THE reference for CAP-Gly recognition of -EEY-COO-; chain A domain, chain C = EEY 266-268 with free carboxylate; motif GKNDG at 67-71; 0 clashes, 138 atom pairs < 4.5 A |
| 2HQH | p150Glued + CLIP-170 | X-ray | - |  |
| 2HKN | p150Glued CAP-Gly | X-ray | - |  |
| 3RDV | SLAIN2c-CLIP CG1 complex | X-ray | - | CAP-Gly + peptides ending in free Tyr |
| 9YNE | dynein-1 motor + dynactin-p150glued | EM | - | megacomplex, no free peptide |
| 9YNH | full-length dynein-1 + dynactin | EM | - | as above |
| 8PR0 | dynein-A + dynactin-p150glued + IC-LC tower | EM | - | as above |
| 1JFF | ab-tubulin, zinc sheets | electron crystallography | 3.5 A | canonical alpha-tubulin numbering; chain A resolved 2-439, so the C-terminal tail is unresolved |
| 3RYF | GTP-tubulin : RB3 stathmin-like domain | X-ray | 2.52 A |  |
| 5IJ0 | human TUBB3 microtubule | EM | 3.8 A |  |
| 7SJ7 | recombinant human microtubule | EM | 3.8 A |  |

Total: **22 entries**, downloaded as mmCIF.

## EMDB

| entry | description | resolution | notes |
|---|---|---|---|
| EMD-2447 | ternary complex of alpha-tubulin, TBCE and TBCB | 21 A | the ONLY experimental structure of a TBCB-containing complex. 32^3 grid, 4.66 A voxels, ~1.3 voxels per residue. Cannot assign interfaces. |
| EMD-2449 | human TBCE alone | 17 A | not used |

## AlphaFold Database models (v6)

| accession | protein |
|---|---|
| Q9JJZ2 | mouse TUBA8 |
| Q9D1E6 | mouse TBCB |
| Q8CIV8 | mouse TBCE |
| P68369 | mouse TUBA1A |
| P68368 | mouse TUBA4A |
| P05213 | mouse TUBA1B |
| Q99426 | human TBCB |
| Q9NY65 | human TUBA8 |

## AlphaFold 3 predictions supplied by the user

| job | chains | models | ipTM (mean of 5) | MSA depth |
|---|---|---|---|---|
| `fold_tbcb_tuba8` | TBCB 244 aa + TUBA8 449 aa | 5 | **0.8300** | 50 000 paired |
| `fold_tbcb_tuba1a` | TBCB 244 aa + TUBA1A 451 aa | 5 | **0.8160** | 50 000 paired |
| `fold_tuba1a` | TUBA1A 451 aa | 5 | n/a (monomer, pTM 0.92) | 50 000 paired |
| `fold_tuba1b` | TUBA1B 451 aa | 5 | n/a (monomer, pTM 0.93) | 50 000 paired |
| `fold_tuba1a_31_50` | 20-mer QPDGQMPSDKTIGGGDDSFN | 5 | n/a (pTM 0.07, 75% disordered) | - |

Each job supplied full PAE matrices, per-model confidences, templates and MSAs. The paired
alpha-tubulin MSA (44 191 usable match-state sequences after gap filtering) is the input to the
direct coupling analysis.

## Non-computational input

Proteomics, supplied by the user: TUBA8 (Q9JJZ2) has valid quantification in all 60 tissue TMT
channels, covering cerebellum, hippocampus, prefrontal cortex, striatum and spinal cord. TUBA8
is therefore **not** absent from the unaffected regions, which excludes an expression-floor
explanation of the regional phenotype and makes a TUBA8-specific rate difference necessary.
