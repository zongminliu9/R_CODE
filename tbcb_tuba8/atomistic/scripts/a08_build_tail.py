"""CAP-Gly : alpha-tail systems taken from the AF3 TBCB:alpha complexes (the pose whose groove geometry AF3 reproduces
from 2E4H/2HL3 to 0.2-0.4 A), for relative binding free energies of the terminal residue (Y<->F) and the K188A
double-mutant cycle, and for tauRAMD residence times.

Complex leg : TBCB residues 134-233 (CAP-Gly domain, ACE/NME capped) + the last 8 residues of the alpha chain (ACE
              capped N-terminus, free C-terminal carboxylate as in the protein).
Free leg    : the same capped 8-residue tail alone in water.
Optional    : --k188a mutates TBCB Lys188 to Ala before protonation (PDBFixer), for the double-mutant cycle.
Usage: a08_build_tail.py <sid> <af3 model.cif> --alpha B --tbcb A [--k188a]
Writes build/<sid>_cplx/ and build/<sid>_free/ (solvated.cif.gz, minimized.cif.gz, build.json, minimize.json)
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
ap.add_argument('sid'); ap.add_argument('src'); ap.add_argument('--alpha', default='B'); ap.add_argument('--tbcb', default='A')
ap.add_argument('--k188a', action='store_true'); ap.add_argument('--seg', default='134-233'); ap.add_argument('--ntail', type=int, default=8)
ap.add_argument('--tail_mut', default=None, help='e.g. TYR-451-PHE, applied to the alpha chain before protonation')
ap.add_argument('--legs', default='cplx,free')
a = ap.parse_args()
ff = A1.forcefield()
st = gemmi.read_structure(a.src); st.setup_entities(); st.remove_hydrogens()
out = gemmi.Structure(); out.add_model(gemmi.Model('1'))
for cid in (a.tbcb, a.alpha):
    nc = gemmi.Chain(cid)
    for r in st[0][cid].get_polymer():
        nc.add_residue(r)
    out[0].add_chain(nc)
out.setup_entities()
fx = PDBFixer(pdbfile=io.StringIO(out.make_pdb_string()))
fx.missingResidues = {}
if a.k188a:
    fx.applyMutations(['LYS-188-ALA'], a.tbcb)
if a.tail_mut:
    fx.applyMutations([a.tail_mut], a.alpha)
fx.findMissingAtoms(); fx.addMissingAtoms()
m = app.Modeller(fx.topology, fx.positions)
m.addHydrogens(ff, pH=7.0)
top, P = m.topology, np.array(m.positions.value_in_unit(C.NM))
lo, hi = [int(x) for x in a.seg.split('-')]
chains = {c.id: c for c in top.chains()}
tail_res = list(chains[a.alpha].residues())[-a.ntail:]
last_alpha = int(tail_res[-1].id)


def extract(segments):
    """segments: list of (chain, first, last, cap_n, cap_c); caps built from flanking residues as in a05.pentapeptide."""
    T = app.Topology(); pos = []
    H = app.element.hydrogen; Ce = app.element.carbon; Oe = app.element.oxygen; Ne = app.element.nitrogen
    def at(r, n): return [x for x in r.atoms() if x.name == n][0]
    def mH(c, partners):
        return [c + 0.109 * (p - c) / np.linalg.norm(p - c) for p in partners]
    for cid, first, last, capn, capc in segments:
        res = {int(r.id): r for r in chains[cid].residues()}
        ch = T.addChain(cid)
        prevC = None
        if capn:
            rA = res[first - 1]; ace = T.addResidue('ACE', ch, str(first - 1))
            ca = P[at(rA, 'CA').index]
            partners = [P[at(rA, 'N').index]] + [P[x.index] for x in rA.atoms() if x.name in ('HA', 'CB', 'HA2', 'HA3')]
            c1 = T.addAtom('C', Ce, ace); pos.append(P[at(rA, 'C').index])
            o1 = T.addAtom('O', Oe, ace); pos.append(P[at(rA, 'O').index])
            cm = T.addAtom('CH3', Ce, ace); pos.append(ca); T.addBond(c1, o1); T.addBond(c1, cm)
            for k, x in enumerate(mH(ca, partners[:3])):
                h = T.addAtom(f'HH3{k + 1}', H, ace); pos.append(x); T.addBond(cm, h)
            prevC = c1
        keep = {}
        for k in range(first, last + 1):
            r = res[k]; rr = T.addResidue(r.name, ch, r.id)
            for x in r.atoms():
                if not capn and k == first and x.name == 'H':
                    continue
                keep[x.index] = T.addAtom(x.name, x.element, rr); pos.append(P[x.index])
            if prevC is not None:
                T.addBond(prevC, keep[at(r, 'N').index])
            prevC = keep[at(r, 'C').index]
        for b1, b2 in top.bonds():
            if b1.index in keep and b2.index in keep:
                T.addBond(keep[b1.index], keep[b2.index])
        if capc:
            rN = res[last + 1]; nme = T.addResidue('NME', ch, str(last + 1))
            n1 = T.addAtom('N', Ne, nme); pos.append(P[at(rN, 'N').index])
            h1 = T.addAtom('H', H, nme); pos.append(P[at(rN, 'H').index])
            ca = P[at(rN, 'CA').index]
            partners = [P[at(rN, 'C').index]] + [P[x.index] for x in rN.atoms() if x.name in ('HA', 'CB', 'HA2', 'HA3')]
            cm = T.addAtom('CH3', Ce, nme); pos.append(ca)
            T.addBond(prevC, n1); T.addBond(n1, h1); T.addBond(n1, cm)
            for k, x in enumerate(mH(ca, partners[:3])):
                h = T.addAtom(f'HH3{k + 1}', H, nme); pos.append(x); T.addBond(cm, h)
    return T, np.array(pos)


tail_seg = (a.alpha, int(tail_res[0].id), last_alpha, True, False)
for leg, segs in (('cplx', [(a.tbcb, lo, hi, True, True), tail_seg]), ('free', [tail_seg])):
    if leg not in a.legs.split(','):
        continue
    W = os.path.join(A1.ROOT, 'build', f'{a.sid}_{leg}'); os.makedirs(W, exist_ok=True)
    T, X = extract(segs)
    mm_ = app.Modeller(T, X * C.NM)
    ext = float(np.max(X.max(0) - X.min(0)))
    if leg == 'cplx':
        mm_.addSolvent(ff, model='tip3p', padding=1.2 * C.NM, boxShape='dodecahedron', ionicStrength=0.15 * unit.molar, neutralize=True)
    else:
        L = max(ext + 2.4, 4.0)
        mm_.addSolvent(ff, model='tip3p', boxSize=mm.Vec3(L, L, L) * C.NM, ionicStrength=0.15 * unit.molar, neutralize=True)
    with gzip.open(os.path.join(W, 'solvated.cif.gz'), 'wt') as fh:
        app.PDBxFile.writeFile(mm_.topology, mm_.positions, fh, keepIds=True)
    seq = {c.id: ''.join(A1.AA3.get(r.name, r.name) and {v: k for k, v in A1.AA3.items()}.get(r.name, 'x') for r in c.residues()) for c in mm_.topology.chains() if c.id in (a.alpha, a.tbcb)}
    meta = dict(sid=f'{a.sid}_{leg}', src=os.path.basename(a.src), k188a=a.k188a, tail_mut=a.tail_mut, tbcb_segment=a.seg if leg == 'cplx' else None,
                tail=[int(tail_res[0].id), last_alpha], sequences=seq, n_atoms=mm_.topology.getNumAtoms())
    json.dump(meta, open(os.path.join(W, 'build.json'), 'w'), indent=1)
    print(json.dumps(meta))
    extra = ['--chains-a', a.tbcb, '--chains-b', a.alpha] if leg == 'cplx' else []
    subprocess.run([sys.executable, '-W', 'ignore', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a02_minimize.py'), W] + extra, check=True)
