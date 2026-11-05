"""TBCE CAP-Gly : alpha-tail complex for the transfer-selectivity cycle.

No structure has an alpha tail in the TBCE groove (AF3 leaves it empty; 9M1L/9M1M resolve alpha only to 437/438),
so the tail pose is transferred from the TBCB groove of the same AF3 model: the TBCE CAP-Gly of 9M1L chain E is
superposed on TBCB 134-233 using sequence-aligned CA pairs (gemmi global alignment, BLOSUM62), with iterative
rejection of pairs farther than 2 A, and the tail coordinates are copied across.  The groove residues equivalent to
TBCB F187/K188/K203/N204/F223 are read from the same alignment and stored for the pose restraint.
Usage: a08b_build_tbce_tail.py <sid> <af3 TBCB:alpha model.cif> [--tbce-range 1-110]
"""
import os, sys, io, json, gzip, argparse, subprocess
import numpy as np
import gemmi
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a01_build as A1
import common as C

ap = argparse.ArgumentParser()
ap.add_argument('sid'); ap.add_argument('src'); ap.add_argument('--tbce-range', default='1-110')
ap.add_argument('--ntail', type=int, default=8)
a = ap.parse_args()
ff = A1.forcefield()
af = gemmi.read_structure(a.src); af.remove_hydrogens()
tb = [r for r in af[0]['A'].get_polymer() if 134 <= r.seqid.num <= 233]
tail = list(af[0]['B'].get_polymer())[-(a.ntail + 1):]          # one extra residue for the ACE cap
ex = gemmi.read_structure(os.path.join(A1.ROOT, 'inputs', '9M1L.cif')); ex.remove_hydrogens()
lo, hi = [int(x) for x in a.tbce_range.split('-')]
te = [r for r in ex[0]['E'].get_polymer() if lo <= r.seqid.num <= hi]
s1 = gemmi.one_letter_code([r.name for r in tb]); s2 = gemmi.one_letter_code([r.name for r in te])
from Bio import Align
from Bio.Align import substitution_matrices
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -10.0; aligner.extend_gap_score = -0.5
aligner.end_gap_score = 0.0; aligner.mode = 'global'
aln_ = aligner.align(s1, s2)[0]
pairs = []
for (a0, a1), (b0, b1) in zip(*aln_.aligned):
    for k in range(a1 - a0):
        pairs.append((a0 + k, b0 + k))
ident = sum(1 for p, q in pairs if s1[p] == s2[q])
def ca(r): return np.array(r['CA'][0].pos.tolist())
X = np.array([ca(te[q]) for p, q in pairs]); Y = np.array([ca(tb[p]) for p, q in pairs])
keep = np.ones(len(pairs), bool)
for it in range(6):
    xc, yc = X[keep].mean(0), Y[keep].mean(0)
    U, S, Vt = np.linalg.svd((X[keep] - xc).T @ (Y[keep] - yc)); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Z = (R @ (X - xc).T).T + yc
    dev = np.linalg.norm(Z - Y, axis=1)
    keep = dev < 2.0
rmsd = float(np.sqrt(np.mean(dev[keep] ** 2)))
move = lambda x: (R @ (np.asarray(x) - xc)) + yc       # TBCE frame -> TBCB frame (the tail stays put)
# structural equivalents: TBCE residue whose CA (after superposition) is nearest to each TBCB groove CA
te_ca = np.array([move(ca(r)) for r in te])
struct_equiv = {}
for r in tb:
    if r.seqid.num in (187, 188, 203, 204, 222, 223):
        dd = np.linalg.norm(te_ca - ca(r), axis=1); k = int(np.argmin(dd))
        struct_equiv[r.seqid.num] = (te[k].seqid.num, te[k].name, round(float(dd[k]), 2))
equiv = {}
for p, q in pairs:
    if tb[p].seqid.num in (187, 188, 203, 204, 222, 223):
        equiv[tb[p].seqid.num] = (te[q].seqid.num, te[q].name)
# assemble TBCE segment (moved into the TBCB frame) + tail, heavy atoms
out = gemmi.Structure(); out.add_model(gemmi.Model('1'))
ce = gemmi.Chain('A')
for r in te:
    nr = gemmi.Residue(); nr.name = r.name; nr.seqid = r.seqid; nr.het_flag = 'A'
    for at in r:
        na = gemmi.Atom(); na.name = at.name; na.element = at.element; na.b_iso = 20.0
        na.pos = gemmi.Position(*move(at.pos.tolist())); nr.add_atom(na)
    ce.add_residue(nr)
out[0].add_chain(ce)
cb = gemmi.Chain('B')
for r in tail:
    cb.add_residue(r)
out[0].add_chain(cb)
out.setup_entities()
fx = PDBFixer(pdbfile=io.StringIO(out.make_pdb_string()))
fx.missingResidues = {}
fx.findMissingAtoms(); fx.addMissingAtoms()
m = app.Modeller(fx.topology, fx.positions)
m.addHydrogens(ff, pH=7.0)
# reuse the capping extractor of a08 by writing an intermediate and calling it with the same conventions
sys.argv = ['x']
import importlib.util
spec = importlib.util.spec_from_file_location('a08', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a08_build_tail.py'))
top, P = m.topology, np.array(m.positions.value_in_unit(C.NM))
chains = {c.id: c for c in top.chains()}
resA = {int(r.id): r for r in chains['A'].residues()}
firstA, lastA = min(resA) + 1, max(resA) - 1
resB = [int(r.id) for r in chains['B'].residues()]
from a05_fep import pentapeptide  # noqa (keeps the cap geometry identical to the other builds)
def at(r, n): return [x for x in r.atoms() if x.name == n][0]
T = app.Topology(); pos = []
Hh = app.element.hydrogen; Ce = app.element.carbon; Oe = app.element.oxygen; Ne = app.element.nitrogen
def mH(c, partners): return [c + 0.109 * (p - c) / np.linalg.norm(p - c) for p in partners]
def seg(cid, res, first, last, capn, capc):
    ch = T.addChain(cid); prevC = None
    if capn:
        rA = res[first - 1]; ace = T.addResidue('ACE', ch, str(first - 1)); ca_ = P[at(rA, 'CA').index]
        partners = [P[at(rA, 'N').index]] + [P[x.index] for x in rA.atoms() if x.name in ('HA', 'CB', 'HA2', 'HA3')]
        c1 = T.addAtom('C', Ce, ace); pos.append(P[at(rA, 'C').index]); o1 = T.addAtom('O', Oe, ace); pos.append(P[at(rA, 'O').index])
        cm = T.addAtom('CH3', Ce, ace); pos.append(ca_); T.addBond(c1, o1); T.addBond(c1, cm)
        for k, x in enumerate(mH(ca_, partners[:3])):
            h = T.addAtom(f'HH3{k + 1}', Hh, ace); pos.append(x); T.addBond(cm, h)
        prevC = c1
    kp = {}
    for k in range(first, last + 1):
        r = res[k]; rr = T.addResidue(r.name, ch, r.id)
        for x in r.atoms():
            kp[x.index] = T.addAtom(x.name, x.element, rr); pos.append(P[x.index])
        if prevC is not None:
            T.addBond(prevC, kp[at(r, 'N').index])
        prevC = kp[at(r, 'C').index]
    for b1, b2 in top.bonds():
        if b1.index in kp and b2.index in kp:
            T.addBond(kp[b1.index], kp[b2.index])
    if capc:
        rN = res[last + 1]; nme = T.addResidue('NME', ch, str(last + 1))
        n1 = T.addAtom('N', Ne, nme); pos.append(P[at(rN, 'N').index]); h1 = T.addAtom('H', Hh, nme); pos.append(P[at(rN, 'H').index])
        ca_ = P[at(rN, 'CA').index]
        partners = [P[at(rN, 'C').index]] + [P[x.index] for x in rN.atoms() if x.name in ('HA', 'CB', 'HA2', 'HA3')]
        cm = T.addAtom('CH3', Ce, nme); pos.append(ca_); T.addBond(prevC, n1); T.addBond(n1, h1); T.addBond(n1, cm)
        for k, x in enumerate(mH(ca_, partners[:3])):
            h = T.addAtom(f'HH3{k + 1}', Hh, nme); pos.append(x); T.addBond(cm, h)
seg('A', resA, firstA, lastA, True, True)
resBd = {int(r.id): r for r in chains['B'].residues()}
seg('B', resBd, resB[1], resB[-1], True, False)
X = np.array(pos)
tailI = [i for i, x in enumerate(T.atoms()) if x.residue.chain.id == 'B' and x.element.symbol != 'H']
recI = [i for i, x in enumerate(T.atoms()) if x.residue.chain.id == 'A' and x.element.symbol != 'H']
dmin = float(np.min(np.linalg.norm(X[tailI][:, None] - X[recI][None], axis=2))) * 10
W = os.path.join(A1.ROOT, 'build', f'{a.sid}_cplx'); os.makedirs(W, exist_ok=True)
mm_ = app.Modeller(T, X * C.NM)
mm_.addSolvent(ff, model='tip3p', padding=1.2 * C.NM, boxShape='dodecahedron', ionicStrength=0.15 * unit.molar, neutralize=True)
with gzip.open(os.path.join(W, 'solvated.cif.gz'), 'wt') as fh:
    app.PDBxFile.writeFile(mm_.topology, mm_.positions, fh, keepIds=True)
meta = dict(sid=f'{a.sid}_cplx', src=os.path.basename(a.src), tbce='9M1L chain E', tbce_range=a.tbce_range,
            superposition_ca_rmsd_A=round(rmsd, 3), n_pairs_kept=int(keep.sum()), n_aligned=len(pairs), n_identical=ident,
            alignment=str(aln_),
            tbcb_to_tbce_equivalents_by_alignment=equiv, tbcb_to_tbce_equivalents_by_structure=struct_equiv, min_tail_receptor_heavy_dist_A=round(dmin, 2), tail=[resB[1], resB[-1]],
            n_atoms=mm_.topology.getNumAtoms())
json.dump(meta, open(os.path.join(W, 'build.json'), 'w'), indent=1)
print(json.dumps(meta))
subprocess.run([sys.executable, '-W', 'ignore', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a02_minimize.py'), W,
                '--chains-a', 'A', '--chains-b', 'B'], check=True)
