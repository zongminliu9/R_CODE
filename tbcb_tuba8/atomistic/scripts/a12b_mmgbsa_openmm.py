"""Single-trajectory MM/GBSA on the GPU next to the trajectories (Sherlock): OBC2 (igb=5) with mbondi2 radii from the
ParmEd topologies written by a12_mmpbsa.py prep, 0.15 M salt (Debye-Hueckel kappa), ACE nonpolar term, no cutoff.
dG_GB(frame) = E(complex) - E(receptor) - E(ligand) at the complex geometry (internal terms cancel exactly).
Reported per replica (mean over frames after start_ns) and across replicas (mean, SD): a ranking score, not a Kd.
Usage: a12b_mmgbsa_openmm.py <mmpbsa_dir> <md_system_dir> <out.json> [--start_ns 10] [--every 5]
"""
import os, sys, json, glob, argparse
import numpy as np
import openmm as mm
from openmm import app, unit
import parmed
import mdtraj as md

ap = argparse.ArgumentParser(); ap.add_argument('mm'); ap.add_argument('mdsys'); ap.add_argument('out')
ap.add_argument('--start_ns', type=float, default=10.0); ap.add_argument('--every', type=int, default=5)
ap.add_argument('--dt_ps', type=float, default=20.0); ap.add_argument('--platform', default='CUDA')
a = ap.parse_args()
info = json.load(open(os.path.join(a.mm, 'prep.json')))
cpx = parmed.load_file(os.path.join(a.mm, 'complex.prmtop'))
def mask_atoms(mask):
    sel = parmed.amber.AmberMask(cpx, mask).Selection()
    return np.where(np.array(sel) > 0)[0]
ir = mask_atoms(info['receptor_mask']); il = mask_atoms(info['ligand_mask'])
assert len(ir) + len(il) == len(cpx.atoms)
def make(prm):
    s = parmed.load_file(prm)
    sysm = s.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.OBC2,
                          implicitSolventSaltConc=0.15 * unit.molar)
    plat = mm.Platform.getPlatformByName(a.platform)
    return mm.Context(sysm, mm.VerletIntegrator(0.001), plat, {'Precision': 'mixed'} if a.platform in ('CUDA', 'OpenCL') else {})
# Sherlock GPUs refuse a second CUDA context in one process ("No compatible CUDA device", 2026-11-07), so each of the
# three systems gets its own context in turn, which is released before the next one is made.
SHIFTS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
def mic(d, box):
    """triclinic minimum image of one vector (box rows a, b, c in OpenMM reduced form)"""
    d = d - np.round(d[2] / box[2, 2]) * box[2]; d = d - np.round(d[1] / box[1, 1]) * box[1]; d = d - np.round(d[0] / box[0, 0]) * box[0]
    cand = d[None, :] + SHIFTS @ box
    return cand[np.argmin((cand ** 2).sum(1))]
frames = {}
for rdir in sorted(glob.glob(os.path.join(a.mdsys, 'r*'))):
    segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd')))
    if not segs:
        continue
    t = md.join([md.load_dcd(s, top=os.path.join(rdir, 'solute.pdb')) for s in segs], check_topology=False)
    t = t[int(a.start_ns * 1000 / a.dt_ps)::a.every]
    assert t.n_atoms == len(cpx.atoms), (t.n_atoms, len(cpx.atoms))
    # same atoms in the same order, not just the same count (2026-11-07: topologies prepared before the 11-05 rebuilds)
    tp = md.load_pdb(os.path.join(rdir, 'solute.pdb'), standard_names=False).topology
    bad = [(i, (x.residue.name, x.name), (y.residue.name, y.name)) for i, (x, y) in enumerate(zip(tp.atoms, cpx.atoms))
           if (x.residue.name.strip().upper(), x.name.strip().upper()) != (y.residue.name.strip().upper(), y.name.strip().upper())]
    assert not bad, f'{len(bad)} atoms differ between solute.pdb and complex.prmtop, first {bad[:3]}'
    X = t.xyz.astype(float).copy(); n_shift = 0
    if t.unitcell_vectors is not None:
        for f in range(t.n_frames):
            d = X[f, il].mean(0) - X[f, ir].mean(0)
            L = mic(d, t.unitcell_vectors[f]) - d
            if np.linalg.norm(L) > 1e-6:
                X[f, il] += L; n_shift += 1
    frames[os.path.basename(rdir)] = (X, n_shift)
def energies(prm, sel):
    ctx = make(prm)
    out = {}
    for k, (X, _) in frames.items():
        e = []
        for f in range(len(X)):
            ctx.setPositions(X[f][sel] if sel is not None else X[f])
            e.append(ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
        out[k] = np.array(e)
    del ctx
    import gc; gc.collect()
    return out
Ec = energies(os.path.join(a.mm, 'complex.prmtop'), None)
Er = energies(os.path.join(a.mm, 'receptor.prmtop'), ir)
El = energies(os.path.join(a.mm, 'ligand.prmtop'), il)
res = {}
for k in frames:
    vals = Ec[k] - Er[k] - El[k]
    res[k] = dict(n_frames=int(len(vals)), mean=float(vals.mean()), sd_frames=float(vals.std(ddof=1)),
                  sem_naive=float(vals.std(ddof=1) / np.sqrt(len(vals))), frames_ligand_reimaged=int(frames[k][1]))
    print(k, res[k], flush=True)
m = [v['mean'] for v in res.values()]
out = dict(system=os.path.basename(a.mdsys.rstrip('/')), GB=dict(mean=float(np.mean(m)) if m else None,
           sd_replicas=float(np.std(m, ddof=1)) if len(m) > 1 else None, n_replicas=len(m)), replicas=res,
           method='single-trajectory MM/GBSA, OBC2 mbondi2, 0.15 M, no entropy; ranking score, not an affinity')
json.dump(out, open(a.out, 'w'), indent=1)
print(json.dumps(out['GB']))
