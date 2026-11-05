"""Residue-level physics from the unbiased MD: for every TUBA8-divergent position and the 30-50 region of one system,
per replica (mean over frames after 10 ns), then mean +/- SD over replicas:
  relative SASA (Shrake-Rupley, Tien et al. 2013 maximum ASA), backbone+side-chain H-bonds (Baker-Hubbard, occupancy
  >= 0.10), salt bridges (carboxylate O to Lys NZ / Arg NE,NH1,NH2 within 4 A, occupancy), hydrophobic contacts
  (side-chain C-C < 4.5 A to non-neighbouring residues), and the residue's interaction energy with the rest of the
  protein + GTP.Mg (Coulomb, dielectric 1, plus Lennard-Jones, 1.2 nm cutoff, residues |i-j| >= 2), using the
  simulated force-field parameters.  Interaction energies are potential-energy terms at sampled geometries, not
  free energies.
Usage: a13_residue_analysis.py <md_root>/<sid> <build_dir> <out.json> [--skip_ns 10] [--stride 5]
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
import openmm as mm
from openmm import app, unit

ap = argparse.ArgumentParser()
ap.add_argument('mdsys'); ap.add_argument('build'); ap.add_argument('out')
ap.add_argument('--skip_ns', type=float, default=10.0); ap.add_argument('--stride', type=int, default=5)
ap.add_argument('--dt_ps', type=float, default=20.0)
ap.add_argument('--max_ns', type=float, default=None, help='use frames up to this time only (common window for interim comparisons)')
a = ap.parse_args()
def bse(x, nb=4):
    """standard error of the mean from nb contiguous blocks (within one replica)"""
    x = np.asarray(x, float)
    if len(x) < 2 * nb:
        return None
    m = [b.mean() for b in np.array_split(x, nb)]
    return float(np.std(m, ddof=1) / np.sqrt(nb))


MAXASA = {'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP': 193, 'CYS': 167, 'GLN': 225, 'GLU': 223, 'GLY': 104, 'HIS': 224,
          'HID': 224, 'HIE': 224, 'HIP': 224, 'ILE': 197, 'LEU': 201, 'LYS': 236, 'MET': 224, 'PHE': 240, 'PRO': 159,
          'SER': 155, 'THR': 172, 'TRP': 285, 'TYR': 263, 'VAL': 174}
HYD = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'TYR', 'CYS'}
DIVERGENT = [7, 16, 24, 32, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 50, 54, 58, 67, 75, 78, 80, 82, 110, 114, 117,
             126, 128, 130, 137, 159, 170, 232, 270, 275, 293, 294, 299, 300, 317, 330, 334, 340, 425, 437]
POS = sorted(set(DIVERGENT) | set(range(30, 51)))
# force-field parameters for the solute (same order as the saved DCD selection)
top_full, pos_full = C.load_cif(os.path.join(a.build, 'minimized.cif.gz'))
ff = C.forcefield()
sysf = ff.createSystem(top_full, nonbondedMethod=app.NoCutoff, constraints=None)
nb = [f for f in sysf.getForces() if isinstance(f, mm.NonbondedForce)][0]
par = np.array([[nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge),
                 nb.getParticleParameters(i)[1].value_in_unit(unit.nanometer),
                 nb.getParticleParameters(i)[2].value_in_unit(unit.kilojoule_per_mole)] for i in range(sysf.getNumParticles())])
reps = sorted(glob.glob(os.path.join(a.mdsys, 'r*')))
per_rep = []
windows = []
for rdir in reps:
    segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd')))
    pdb = os.path.join(rdir, 'solute.pdb')
    if not segs:
        continue
    t = md.join([md.load_dcd(s, top=pdb) for s in segs], check_topology=False)
    t_end_ns = t.n_frames * a.dt_ps / 1000.0
    t = t[int(a.skip_ns * 1000 / a.dt_ps):(int(a.max_ns * 1000 / a.dt_ps) if a.max_ns else None):a.stride]
    if t.n_frames < 10:
        continue
    top = t.topology
    # map solute atoms of the DCD to full-system indices (same order as a03_md.py save_idx)
    shell = C.mg_shell_pairs(top_full, pos_full)
    atoms_full = list(top_full.atoms())
    save_idx = sorted(set(C.solute_indices(top_full)) | {b for _, b in shell} |
                      {x.index for _, b in shell for x in atoms_full[b].residue.atoms()})
    assert len(save_idx) == top.n_atoms, (len(save_idx), top.n_atoms)
    P = par[save_idx]
    res = {r.resSeq: r for r in top.residues if r.is_protein}
    sasa = md.shrake_rupley(t, mode='residue')
    hb = md.baker_hubbard(t, freq=0.10, periodic=False)
    hb_frames = None
    out = {}
    xyz = t.xyz
    prot_res = [r for r in top.residues if r.is_protein]
    other = [r for r in top.residues]
    for p in POS:
        if p not in res:
            continue
        r = res[p]; ri = [x.index for x in r.atoms]
        # md.shrake_rupley returns nm^2; MAXASA is in A^2 (x100), and the result is a percentage (x100)
        rel = sasa[:, r.index] * 1.0e4 / MAXASA.get(r.name, 200)
        rec = {'resname': r.name, 'rel_sasa_mean': float(np.mean(rel)), 'rel_sasa_block_se': bse(rel),
               'sasa_A2_mean': float(np.mean(sasa[:, r.index]) * 100.0)}
        rec['hbonds'] = [f'{top.atom(d)}-{top.atom(acc)}' for d, h, acc in hb if top.atom(d).residue.index == r.index or top.atom(acc).residue.index == r.index]
        # salt bridges
        sb = []
        acid = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}; base = {'LYS': ('NZ',), 'ARG': ('NE', 'NH1', 'NH2')}
        def gi(rr, names): return [x.index for x in rr.atoms if x.name in names]
        mine = gi(r, acid.get(r.name, ()) + base.get(r.name, ()))
        if mine:
            partners = [rr for rr in prot_res if rr.index != r.index and ((r.name in acid and rr.name in base) or (r.name in base and rr.name in acid))]
            for rr in partners:
                th = gi(rr, acid.get(rr.name, ()) + base.get(rr.name, ()))
                pairs = np.array([[i, j] for i in mine for j in th])
                d = md.compute_distances(t, pairs, periodic=False).min(axis=1)
                occ = float(np.mean(d < 0.40))
                if occ >= 0.10:
                    sb.append((f'{rr.name}{rr.resSeq}', round(occ, 3)))
        rec['salt_bridges'] = sb
        # hydrophobic contacts
        if r.name in HYD:
            scC = [x.index for x in r.atoms if x.element.symbol == 'C' and x.name not in ('C', 'CA')]
            hc = 0.0
            for rr in prot_res:
                if abs(rr.resSeq - p) < 2 or rr.name not in HYD:
                    continue
                oc = [x.index for x in rr.atoms if x.element.symbol == 'C' and x.name not in ('C', 'CA')]
                if not oc or not scC:
                    continue
                d = md.compute_distances(t[::4], np.array([[i, j] for i in scC for j in oc]), periodic=False)
                hc += float(np.mean((d < 0.45).any(axis=1)))
            rec['hydrophobic_partners_mean'] = hc
        # interaction energy with the rest (|i-j|>=2 residues, GTP, Mg, shell waters)
        rest = [x.index for rr in other for x in rr.atoms if not (rr.is_protein and abs(rr.resSeq - p) < 2)]
        qi, si, ei = P[ri, 0], P[ri, 1], P[ri, 2]
        qj, sj, ej = P[rest, 0], P[rest, 1], P[rest, 2]
        Ec, El = [], []
        for f in range(0, t.n_frames, 2):
            d = np.linalg.norm(xyz[f, ri][:, None] - xyz[f, rest][None], axis=2)
            m = d < 1.2
            sig = 0.5 * (si[:, None] + sj[None]); eps = np.sqrt(ei[:, None] * ej[None])
            Ec.append(np.sum(np.where(m, 138.935458 * qi[:, None] * qj[None] / d, 0.0)))
            sr = (sig / d) ** 6
            El.append(np.sum(np.where(m, 4 * eps * (sr * sr - sr), 0.0)))
        rec['E_coul_kcal'] = float(np.mean(Ec)) / 4.184; rec['E_lj_kcal'] = float(np.mean(El)) / 4.184
        se_c, se_l = bse(Ec), bse(El)
        rec['E_coul_block_se_kcal'] = se_c / 4.184 if se_c is not None else None
        rec['E_lj_block_se_kcal'] = se_l / 4.184 if se_l is not None else None
        out[p] = rec
    per_rep.append(out)
    windows.append([a.skip_ns, min(a.max_ns, t_end_ns) if a.max_ns else t_end_ns, int(t.n_frames)])
    print(os.path.basename(rdir), t.n_frames, 'frames', flush=True)
summ = {}
for p in POS:
    vals = [o[p] for o in per_rep if p in o]
    if not vals:
        continue
    s = {'resname': vals[0]['resname'], 'n_rep': len(vals)}
    for k in ('rel_sasa_mean', 'E_coul_kcal', 'E_lj_kcal', 'hydrophobic_partners_mean'):
        v = [x[k] for x in vals if k in x]
        if v:
            s[k] = [float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else None]
    s['block_se_per_rep'] = [{k: x.get(k) for k in ('rel_sasa_block_se', 'E_coul_block_se_kcal', 'E_lj_block_se_kcal')} for x in vals]
    s['salt_bridges'] = [x['salt_bridges'] for x in vals]
    s['hbonds'] = [x['hbonds'] for x in vals]
    summ[p] = s
json.dump({'system': os.path.basename(a.mdsys.rstrip('/')), 'replicas': len(per_rep), 'window_ns_frames_per_rep': windows,
           'per_residue': summ}, open(a.out, 'w'), indent=1)
print('wrote', a.out)
