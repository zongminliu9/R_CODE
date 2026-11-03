#!/usr/bin/env python3
"""Global (Needleman-Wunsch, affine gap, BLOSUM62) alignment of alpha-tubulin isoforms.
Stdlib only -- alpha-tubulins are >95% identical so exact DP global alignment is unambiguous.
Emits: per-pair divergence tables, TUBA8-specific position list, canonical-numbering map."""
import json, os, itertools

B62_RAW = """A R N D C Q E G H I L K M F P S T W Y V
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
lines = B62_RAW.strip().split("\n"); cols = lines[0].split()
B62 = {}
for ln in lines[1:]:
    p = ln.split(); r = p[0]
    for c, v in zip(cols, p[1:]): B62[(r, c)] = int(v)
def sc(a, b): return B62.get((a, b), B62.get((b, a), -4))

def nw_affine(s1, s2, gap_open=-11, gap_ext=-1):
    n, m = len(s1), len(s2); NEG = -10**9
    M = [[NEG]*(m+1) for _ in range(n+1)]; Ix = [[NEG]*(m+1) for _ in range(n+1)]; Iy = [[NEG]*(m+1) for _ in range(n+1)]
    pM = [[0]*(m+1) for _ in range(n+1)]; pIx = [[0]*(m+1) for _ in range(n+1)]; pIy = [[0]*(m+1) for _ in range(n+1)]
    M[0][0] = 0
    for i in range(1, n+1): Ix[i][0] = gap_open + (i-1)*gap_ext; pIx[i][0] = 1
    for j in range(1, m+1): Iy[0][j] = gap_open + (j-1)*gap_ext; pIy[0][j] = 2
    for i in range(1, n+1):
        s1i = s1[i-1]
        for j in range(1, m+1):
            d = sc(s1i, s2[j-1])
            cands = (M[i-1][j-1]+d, Ix[i-1][j-1]+d, Iy[i-1][j-1]+d)
            k = max(range(3), key=lambda t: cands[t]); M[i][j] = cands[k]; pM[i][j] = k
            cands = (M[i-1][j]+gap_open, Ix[i-1][j]+gap_ext)
            k = max(range(2), key=lambda t: cands[t]); Ix[i][j] = cands[k]; pIx[i][j] = (0 if k == 0 else 1)
            cands = (M[i][j-1]+gap_open, Iy[i][j-1]+gap_ext)
            k = max(range(2), key=lambda t: cands[t]); Iy[i][j] = cands[k]; pIy[i][j] = (0 if k == 0 else 2)
    i, j = n, m; state = max(range(3), key=lambda t: (M[n][m], Ix[n][m], Iy[n][m])[t])
    a1 = []; a2 = []
    while i > 0 or j > 0:
        if state == 0:
            a1.append(s1[i-1]); a2.append(s2[j-1]); nxt = pM[i][j]; i -= 1; j -= 1; state = nxt
        elif state == 1:
            a1.append(s1[i-1]); a2.append("-"); nxt = pIx[i][j]; i -= 1; state = nxt
        else:
            a1.append("-"); a2.append(s2[j-1]); nxt = pIy[i][j]; j -= 1; state = nxt
        if i == 0 and j > 0: state = 2
        elif j == 0 and i > 0: state = 1
    return "".join(reversed(a1)), "".join(reversed(a2))

idx = json.load(open("data/meta/sequence_index.json"))
ALPHAS = ["Tuba8","Tuba1a","Tuba1b","Tuba1c","Tuba3a","Tuba4a"]
seqs = {g: idx[f"mouse_{g}"]["sequence"] for g in ALPHAS}
REF = "Tuba1a"   # canonical alpha-tubulin numbering reference (451 aa, == 1JFF/bovine numbering)

os.makedirs("results/alignment", exist_ok=True)
report = {}
print(f"{'pair':22s} {'ident':>7s} {'%id':>6s} {'ndiff':>6s} {'gaps':>5s}")
for a, b in itertools.combinations(ALPHAS, 2):
    A, Bb = nw_affine(seqs[a], seqs[b])
    ident = sum(1 for x, y in zip(A, Bb) if x == y and x != "-")
    ncol = sum(1 for x, y in zip(A, Bb) if x != "-" and y != "-")
    gaps = sum(1 for x, y in zip(A, Bb) if x == "-" or y == "-")
    ndiff = ncol - ident
    report[f"{a}|{b}"] = {"identities": ident, "aligned_cols": ncol, "pct_id": round(100*ident/ncol, 2),
                          "n_substitutions": ndiff, "n_gap_cols": gaps}
    print(f"{a+' vs '+b:22s} {ident:7d} {100*ident/ncol:6.2f} {ndiff:6d} {gaps:5d}")
json.dump(report, open("results/alignment/pairwise_identity_mouse.json","w"), indent=1)

# ---------------- divergent-position enumeration in canonical (TUBA1A) numbering ----------------
def aln_map(ref, qry):
    """Return list of (ref_pos|None, qry_pos|None, ref_aa, qry_aa) over alignment columns."""
    A, B = nw_affine(ref, qry); i = j = 0; out = []
    for x, y in zip(A, B):
        ri = qi = None
        if x != "-": i += 1; ri = i
        if y != "-": j += 1; qi = j
        out.append((ri, qi, x, y))
    return out

# structural region assignment for alpha-tubulin (canonical TUBA1A/1JFF numbering)
REGIONS = [
 (1,  4,  "N-term"),            (5,  9,  "S1"),              (10, 24, "H1"),
 (25, 38, "H1-S2 loop (N-side)"),(39,46,"K40 loop / acetyl loop"),
 (47, 62, "H1-S2 loop (C-side)"),(63, 73, "S2"),             (74, 85, "H2 / T2"),
 (86,100, "S3 / T3"),           (101,110,"H3"),              (111,130,"H3-S4 / T3 loop"),
 (131,143,"S4"),                (144,152,"T4 / loop"),       (153,172,"H4"),
 (173,180,"S5"),                (181,190,"T5 loop"),         (191,205,"H5"),
 (206,221,"S6 / H6"),           (222,243,"H6-H7 loop / H6"), (244,272,"H7 (core helix)"),
 (273,285,"H7-H8 loop / T7"),   (286,295,"H8"),              (296,310,"S7"),
 (311,325,"H9 / S8 region"),    (326,340,"S8-H9 loop (M-loop side)"),
 (341,352,"H9-S9"),             (353,370,"S9 / M-loop"),     (371,381,"S10"),
 (382,396,"H10"),               (397,410,"H10-S10 / S10"),   (411,424,"H11"),
 (425,437,"H11-H12 loop"),      (438,451,"H12 / C-term tail"),
]
def region_of(p):
    for lo, hi, nm in REGIONS:
        if lo <= p <= hi: return nm
    return "beyond-H12 (tail)"

ref = seqs[REF]
rows = []
for g in ALPHAS:
    if g == REF: continue
    for ri, qi, x, y in aln_map(ref, seqs[g]):
        if x != y:
            rows.append({"isoform": g, "canon_pos": ri, "iso_pos": qi,
                         "ref_aa": x, "iso_aa": y, "region": region_of(ri) if ri else "insertion"})

# Which positions are TUBA8-specific (TUBA8 differs from REF *and* from every other isoform)?
by_pos = {}
for g in ALPHAS:
    m = {ri: y for ri, qi, x, y in aln_map(ref, seqs[g]) if ri}
    by_pos[g] = m
canon_len = len(ref)
tuba8_rows = []
for p in range(1, canon_len+1):
    col = {g: by_pos[g].get(p, "-") for g in ALPHAS}
    r = ref[p-1]
    if col["Tuba8"] != r:
        others = [col[g] for g in ALPHAS if g != "Tuba8"]
        uniq = all(col["Tuba8"] != o for o in others)
        tuba8_rows.append({"canon_pos": p, "region": region_of(p), "TUBA1A": r,
                           "TUBA8": col["Tuba8"], "TUBA1B": col["Tuba1b"], "TUBA1C": col["Tuba1c"],
                           "TUBA3A": col["Tuba3a"], "TUBA4A": col["Tuba4a"],
                           "tuba8_unique": uniq})

import csv
with open("results/alignment/tuba8_divergent_positions.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(tuba8_rows[0].keys())); w.writeheader(); w.writerows(tuba8_rows)
with open("results/alignment/all_isoform_substitutions.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

nuq = sum(1 for r in tuba8_rows if r["tuba8_unique"])
print(f"\nTUBA8 differs from TUBA1A at {len(tuba8_rows)} canonical positions; {nuq} are TUBA8-UNIQUE "
      f"(no other mouse alpha-tubulin shares them).\n")
hdr = f"{'pos':>4s} {'1A':>2s}>{'8':<2s} {'uniq':>4s}  {'1B':>2s} {'1C':>2s} {'3A':>2s} {'4A':>2s}  region"
print(hdr); print("-"*len(hdr)+"---------")
for r in tuba8_rows:
    print(f"{r['canon_pos']:4d} {r['TUBA1A']:>2s}>{r['TUBA8']:<2s} {'***' if r['tuba8_unique'] else '':>4s}  "
          f"{r['TUBA1B']:>2s} {r['TUBA1C']:>2s} {r['TUBA3A']:>2s} {r['TUBA4A']:>2s}  {r['region']}")

# ---- verify the specific claims from the prior analysis ----
print("\n=== verification of pre-named residues ===")
a8 = seqs["Tuba8"]; a1a = seqs["Tuba1a"]
m8 = {ri: (x, y) for ri, qi, x, y in aln_map(a1a, a8) if ri}
for p, claim in [(40, "A40 (K40 acetylation site)"), (110, "I110V"), (270, "A270V")]:
    x, y = m8.get(p, ("?", "?"))
    print(f"  canon {p:3d}: TUBA1A={x} TUBA8={y}   claim='{claim}'  -> {'CONSISTENT' if y!=x else 'IDENTICAL (claim wrong)'}")
print(f"  C-termini:  TUBA1A ...{a1a[-8:]}   TUBA8 ...{a8[-8:]}   TUBA4A ...{seqs['Tuba4a'][-8:]}")
json.dump(tuba8_rows, open("results/alignment/tuba8_divergent_positions.json","w"), indent=1)
