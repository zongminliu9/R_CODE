"""Repair the builds whose protein chirality is wrong (found 2026-11-05 with a29_chirality_check.py).

1. Mirror images from a07_build_complex.py: the solute was turned to its principal axes with the eigenvector matrix
   of np.linalg.eigh, whose determinant can be -1, so the "rotation" was a reflection for c_tbcb_a1a_m0, c_tbcb_a8_m0,
   x_9M1L_tbce_a1a and x_9M1L_tbce_a8 (every protein residue D, GTP mirrored as well).  The solvated and minimised
   structures are reflected back (x -> Lx - x in their rectangular boxes).  This gives the intended L structure with
   the same solvent configuration.  The force field used here (ff14SB from ISOLDE + Meagher GTP + TIP3P / Li-Merz)
   has no CMAP and every torsion phase it applies to these systems is 0 or pi, so the potential energy is invariant
   under reflection: the reflected minimised structure is a minimum of the same energy function.  This is checked
   below by evaluating the potential energy of both coordinate sets before anything is written.
2. t_a1aYF_K188A_cplx: its solvated structure is L, but Phe451 (the C-terminal residue of the tail) became D during
   minimisation.  It is re-minimised with a02_minimize.py, which now also restrains CB in stage 1 and stops when a
   residue changes handedness.
Old files are kept with the suffix _mirror (or _inverted); build.json records what was done.
Usage: a30_fix_chirality.py <build_root>
"""
import os, sys, json, gzip, shutil, subprocess
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

ROOT = sys.argv[1]
MIRRORED = ['c_tbcb_a1a_m0', 'c_tbcb_a8_m0', 'x_9M1L_tbce_a1a', 'x_9M1L_tbce_a8']
INVERTED_ON_MINIMISATION = {'t_a1aYF_K188A_cplx': ['--chains-a', 'A', '--chains-b', 'B']}


def handedness(top, X):
    """fraction of L residues (N-CA-C-CB dihedral < 0) per protein chain"""
    out = {}
    for ch in top.chains():
        s = []
        for r in ch.residues():
            A = {a.name: a.index for a in r.atoms()}
            if r.name == 'GLY' or not all(k in A for k in ('N', 'CA', 'C', 'CB')):
                continue
            p0, p1, p2, p3 = X[A['N']], X[A['CA']], X[A['C']], X[A['CB']]
            b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
            b1 = b1 / np.linalg.norm(b1); v = b0 - b0 @ b1 * b1; w = b2 - b2 @ b1 * b1
            s.append(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)) < 0)
        if s:
            out[ch.id] = float(np.mean(s))
    return out


def energy(system, X, box):
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
    ctx.setPeriodicBoxVectors(*box)
    ctx.setPositions(X * unit.nanometer)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    del ctx
    return e


ff = C.forcefield()
report = {}
for sid in MIRRORED:
    W = os.path.join(ROOT, sid)
    rec = {}
    for name in ('minimized.cif.gz', 'solvated.cif.gz'):
        path = os.path.join(W, name)
        top, pos = C.load_cif(path)
        X = np.array(pos.value_in_unit(unit.nanometer))
        box = top.getPeriodicBoxVectors()
        B = np.array([[v[i].value_in_unit(unit.nanometer) for i in range(3)] for v in box])
        assert np.max(np.abs(B - np.diag(np.diag(B)))) < 1e-9, f'{sid}/{name}: box is not rectangular'
        h0 = handedness(top, X)
        old = os.path.join(W, name.replace('.cif.gz', '_mirror.cif.gz'))
        if all(f > 0.95 for f in h0.values()) and os.path.exists(old):
            print(sid, name, 'already L (fixed earlier), skipped', flush=True); rec[name] = 'already fixed'; continue
        assert all(f < 0.05 for f in h0.values()), f'{sid}/{name}: expected a mirror image, got fraction L {h0}'
        Xr = X.copy(); Xr[:, 0] = B[0, 0] - X[:, 0]
        h1 = handedness(top, Xr)
        assert all(f > 0.95 for f in h1.values()), f'{sid}/{name}: reflection did not give L chains: {h1}'
        r = dict(fraction_L_before=h0, fraction_L_after=h1)
        if name == 'minimized.cif.gz':
            system = ff.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=None,
                                     rigidWater=False, ewaldErrorTolerance=5e-4)
            e0, e1 = energy(system, X, box), energy(system, Xr, box)
            r.update(potential_kJmol_mirror=e0, potential_kJmol_reflected=e1, difference_kJmol=e1 - e0)
            assert abs(e1 - e0) < 1.0, f'{sid}: potential energy changed under reflection by {e1 - e0:.3f} kJ/mol'
        if not os.path.exists(old):
            shutil.copy2(path, old)
        tmp = path + '.tmp'
        with gzip.open(tmp, 'wt') as fh:
            app.PDBxFile.writeFile(top, Xr * unit.nanometer, fh, keepIds=True)
        os.replace(tmp, path)
        rec[name] = r
        print(sid, name, json.dumps(r), flush=True)
    bj = os.path.join(W, 'build.json')
    info = json.load(open(bj))
    info['chirality_fix_20261105'] = dict(reason='a07 principal-axes transform was a reflection (det -1)', files=rec,
                                          kept_as='minimized_mirror.cif.gz, solvated_mirror.cif.gz')
    json.dump(info, open(bj, 'w'), indent=1)
    report[sid] = rec

for sid, extra in INVERTED_ON_MINIMISATION.items():
    W = os.path.join(ROOT, sid)
    path = os.path.join(W, 'minimized.cif.gz')
    keep = os.path.join(W, 'minimized_inverted.cif.gz')
    top, pos = C.load_cif(path)
    if os.path.exists(keep) and all(f == 1.0 for f in handedness(top, np.array(pos.value_in_unit(unit.nanometer))).values()):
        print(sid, 'already re-minimised with all residues L, skipped', flush=True); report[sid] = 'already fixed'; continue
    if not os.path.exists(keep):
        shutil.copy2(path, keep)
    if os.path.exists(os.path.join(W, 'minimize.json')) and not os.path.exists(os.path.join(W, 'minimize_inverted.json')):
        shutil.copy2(os.path.join(W, 'minimize.json'), os.path.join(W, 'minimize_inverted.json'))
    subprocess.run([sys.executable, '-W', 'ignore', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a02_minimize.py'), W] + extra,
                   check=True)
    top, pos = C.load_cif(path)
    h = handedness(top, np.array(pos.value_in_unit(unit.nanometer)))
    assert all(f == 1.0 for f in h.values()), f'{sid}: still not all L after re-minimisation: {h}'
    report[sid] = dict(fraction_L_after=h)
    bj = os.path.join(W, 'build.json')
    info = json.load(open(bj))
    info['chirality_fix_20261105'] = dict(reason='Phe451 became D during the first minimisation', remin='a02 with CB '
                                          'restrained in stage 1 and the handedness gate', fraction_L_after=h,
                                          kept_as='minimized_inverted.cif.gz, minimize_inverted.json')
    json.dump(info, open(bj, 'w'), indent=1)
    print(sid, json.dumps(report[sid]), flush=True)
json.dump(report, open(os.path.join(ROOT, 'chirality_fix_20261105.json'), 'w'), indent=1)
print('CHIRALITY_FIX_DONE')
