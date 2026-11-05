"""Chirality of every protein chain in built structures or trajectory frames (L-amino acids vs mirror images).

For each non-Gly residue with N, CA, C, CB: the dihedral N-CA-C-CB is about -124 deg for an L-amino acid and
+124 deg for its mirror image (D).  A chain is reported L when >= 95 % of its residues are L, D when >= 95 % are D,
and MIXED otherwise.  Reads .pdb, .cif or .cif.gz (mmCIF atom_site with Cartn_x/y/z, auth or label ids).
Usage: a29_chirality_check.py [--strict] <file> [<file> ...]   (one line per chain; exit status 1 if any chain is
       D/MIXED, or with --strict if any single residue is D; --strict also lists the D residues)
"""
import gzip, math, sys


def read(path):
    op = gzip.open if path.endswith('.gz') else open
    atoms = []
    with op(path, 'rt') as fh:
        if '.pdb' in path:
            for l in fh:
                if l.startswith(('ATOM', 'HETATM')):
                    atoms.append((l[21], int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                                  (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
            return atoms
        cols, inloop = [], False
        for l in fh:
            if l.startswith('_atom_site.'):
                cols.append(l.strip().split('.', 1)[1]); inloop = True; continue
            if inloop and l.startswith(('ATOM', 'HETATM')):
                f = l.split(); ix = {c: i for i, c in enumerate(cols)}
                ch = f[ix['auth_asym_id']] if 'auth_asym_id' in ix else f[ix['label_asym_id']]
                rs = f[ix['auth_seq_id']] if 'auth_seq_id' in ix else f[ix['label_seq_id']]
                atoms.append((ch, int(rs), f[ix['label_comp_id']], f[ix['label_atom_id']].strip('"'),
                              (float(f[ix['Cartn_x']]), float(f[ix['Cartn_y']]), float(f[ix['Cartn_z']]))))
            elif inloop and atoms and l.startswith(('#', 'loop_')):
                break
    return atoms


def dihedral(p0, p1, p2, p3):
    sub = lambda a, b: (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    cross = lambda a, b: (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    b0, b1, b2 = sub(p0, p1), sub(p2, p1), sub(p3, p2)
    n = math.sqrt(dot(b1, b1)); b1 = (b1[0] / n, b1[1] / n, b1[2] / n)
    v = sub(b0, tuple(dot(b0, b1) * x for x in b1)); w = sub(b2, tuple(dot(b2, b1) * x for x in b1))
    return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))


bad = False
strict = '--strict' in sys.argv
for path in [x for x in sys.argv[1:] if x != '--strict']:
    res = {}
    for ch, rs, rn, an, xyz in read(path):
        res.setdefault((ch, rs, rn), {})[an] = xyz
    per = {}
    for (ch, rs, rn), A in res.items():
        if rn in ('GLY', 'HOH', 'WAT') or not all(k in A for k in ('N', 'CA', 'C', 'CB')):
            continue
        per.setdefault(ch, []).append((dihedral(A['N'], A['CA'], A['C'], A['CB']), f'{rn}{rs}'))
    for ch, v in sorted(per.items()):
        fL = sum(x < 0 for x, _ in v) / len(v)
        call = 'L' if fL >= 0.95 else ('D (mirror image)' if fL <= 0.05 else 'MIXED')
        dl = [n for x, n in v if x > 0]
        bad |= call != 'L' or (strict and bool(dl))
        extra = f'\tD residues: {" ".join(dl)}' if strict and dl and call == 'L' else ''
        print(f'{path}\tchain {ch}\t{len(v)} residues\tfraction L {fL:.3f}\t{call}{extra}')
sys.exit(1 if bad else 0)
