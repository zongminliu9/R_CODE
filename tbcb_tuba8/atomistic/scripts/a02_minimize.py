"""Energy minimisation with a full residual-force record, for every starting structure used in this campaign.

Records, before and after minimisation: total potential energy and its per-force decomposition, the maximum atomic
force (with the atom it acts on), the RMS force over all atoms and over solute atoms, and geometry changes
(backbone RMSD to the start, the Mg2+ inner shell, and for complexes the interface contact count and the RMSD of
interface residues).  Minimised as a fully flexible system (no bond constraints, flexible TIP3P), so the residual force is the
true gradient norm.  Two stages: (1) backbone heavy atoms and CB restrained (k = 1000 kJ/mol/nm^2) and Mg2+ ligand
distances held by flat-bottom restraints, so clashes are relieved by side chains and solvent rather than by moving
the fold; (2) everything free, minimised to a 10 kJ/mol/nm force tolerance.
Handedness (added 2026-11-05 after minimisation inverted strained residues in several builds): stops before
minimising if any protein residue is D; in both stages a flat-bottom restraint on every N-CA-C-CB dihedral
(zero within 35 deg of the L value -124 deg, 1000 kJ/mol/rad^2 outside) prevents inversion without acting on
normal geometry; CB is also restrained in stage 1; and minimized.cif.gz is not written if any residue's
dihedral still changes sign.  Reported energies exclude every restraint.

Usage: a02_minimize.py <build_dir> [--platform CPU|CUDA|OpenCL] [--chains-a A --chains-b B]
Writes <build_dir>/minimized.cif.gz and <build_dir>/minimize.json
"""
import os, sys, json, gzip, argparse, time, math
import numpy as np
import openmm as mm
from openmm import app, unit

HERE = os.path.dirname(os.path.abspath(__file__))
FF = os.path.join(HERE, '..', 'ff')
KJ = unit.kilojoule_per_mole
NM = unit.nanometer
SOLVENT = {'HOH', 'NA', 'CL', 'WAT'}


def forcefield():
    return app.ForceField(os.path.join(FF, 'amberff14SB.xml'), os.path.join(FF, 'gtp_meagher.xml'),
                          os.path.join(FF, 'water_ions_tip3p_iod.xml'))


def load(path):
    app.Topology.loadBondDefinitions(os.path.join(FF, 'gtp_bonds.xml'))
    with gzip.open(path, 'rt') as fh:
        pdbx = app.PDBxFile(fh)
    return pdbx.topology, pdbx.positions


def make_system(top, ff, hmr=False, flexible=False):
    kw = dict(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * NM, constraints=None if flexible else app.HBonds,
              rigidWater=not flexible, ewaldErrorTolerance=5e-4)
    if hmr:
        kw['hydrogenMass'] = 3.0 * unit.amu
    s = ff.createSystem(top, **kw)
    for i, f in enumerate(s.getForces()):
        f.setForceGroup(min(i, 31))
    return s


def snapshot(ctx, system, top, solute, label):
    names = {i: type(f).__name__ for i, f in enumerate(system.getForces())}
    out = {}
    tot = ctx.getState(getEnergy=True, getForces=True)
    E = tot.getPotentialEnergy().value_in_unit(KJ)
    F = tot.getForces(asNumpy=True).value_in_unit(KJ / NM)
    mag = np.linalg.norm(F, axis=1)
    # constrained hydrogens and virtual sites still carry a force here; report over all atoms and over solute
    atoms = list(top.atoms())
    imax = int(np.argmax(mag))
    a = atoms[imax]
    out['label'] = label
    out['potential_kJmol'] = E
    out['potential_kcalmol'] = E / 4.184
    out['max_force_kJmol_nm'] = float(mag[imax])
    out['max_force_kcalmol_A'] = float(mag[imax]) / 41.84
    out['max_force_atom'] = f'{a.residue.chain.id}:{a.residue.name}{a.residue.id}:{a.name}'
    out['rms_force_all_kcalmol_A'] = float(np.sqrt(np.mean(mag ** 2))) / 41.84
    out['rms_force_solute_kcalmol_A'] = float(np.sqrt(np.mean(mag[solute] ** 2))) / 41.84
    out['n_atoms_force_gt_1000kJmolnm'] = int(np.sum(mag > 1000.0))
    comp = {}
    for g, n in names.items():
        if g > 31:
            continue
        e = ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy().value_in_unit(KJ)
        comp[n] = comp.get(n, 0.0) + e / 4.184
    out['components_kcalmol'] = comp
    return out


def handedness(top, X):
    """per protein residue (non-Gly, with N, CA, C, CB): True if L (N-CA-C-CB dihedral < 0)"""
    out = {}
    for r in top.residues():
        A = {a.name: a.index for a in r.atoms()}
        if r.name in {'GLY'} | SOLVENT or not all(k in A for k in ('N', 'CA', 'C', 'CB')):
            continue
        p0, p1, p2, p3 = X[A['N']], X[A['CA']], X[A['C']], X[A['CB']]
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1); v = b0 - b0 @ b1 * b1; w = b2 - b2 @ b1 * b1
        out[(r.chain.id, r.id, r.name)] = bool(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)) < 0)
    return out


def mg_shell(top, pos, cutoff=0.26):
    P = np.array(pos.value_in_unit(NM))
    atoms = list(top.atoms())
    mg = [a.index for a in atoms if a.residue.name == 'MG']
    res = []
    for m in mg:
        d = np.linalg.norm(P - P[m], axis=1)
        idx = [i for i in np.argsort(d) if 0 < d[i] < cutoff]
        res.append([(f'{atoms[i].residue.name}{atoms[i].residue.id}:{atoms[i].name}', round(float(d[i]) * 10, 2)) for i in idx])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('wd')
    ap.add_argument('--platform', default='CPU')
    ap.add_argument('--chains-a', default=None, help='comma list of chain ids for partner A (complexes)')
    ap.add_argument('--chains-b', default=None)
    ap.add_argument('--input', default='solvated.cif.gz')
    a = ap.parse_args()
    t0 = time.time()
    top, pos = load(os.path.join(a.wd, a.input))
    ff = forcefield()
    # fully flexible (no constraints, flexible TIP3P) so that the reported residual force is the true gradient;
    # with constraints the raw force on constrained atoms includes the constraint reaction and is not a residual
    system = make_system(top, ff, flexible=True)
    atoms = list(top.atoms())
    solute = np.array([at.index for at in atoms if at.residue.name not in SOLVENT])
    bb = [at.index for at in atoms if at.name in ('N', 'CA', 'C', 'O') and at.residue.name not in SOLVENT | {'GTP', 'MG'}]
    cb = [at.index for at in atoms if at.name == 'CB' and at.residue.name not in SOLVENT | {'GTP', 'MG'}]
    hand0 = handedness(top, np.array(pos.value_in_unit(NM)))
    dres = [f'{k[0]}:{k[2]}{k[1]}' for k, v in hand0.items() if not v]
    if dres:
        raise SystemExit(f'D residues in the input ({len(dres)}): {dres[:20]} - not minimised')
    # stage-1 restraints
    rest = mm.CustomExternalForce('kbb*periodicdistance(x, y, z, x0, y0, z0)^2')
    rest.addGlobalParameter('kbb', 1000.0)
    for p in ('x0', 'y0', 'z0'):
        rest.addPerParticleParameter(p)
    P0 = pos.value_in_unit(NM)
    for i in bb + cb:
        rest.addParticle(i, P0[i])
    rest.setForceGroup(31)
    system.addForce(rest)
    fb = mm.CustomBondForce('kmg*step(r-r0)*(r-r0)^2')
    fb.addGlobalParameter('kmg', 50000.0); fb.addPerBondParameter('r0')
    Parr = np.array(P0)
    mgs = [at.index for at in atoms if at.residue.name == 'MG']
    shell0 = []
    for m in mgs:
        d = np.linalg.norm(Parr - Parr[m], axis=1)
        for i in np.where((d > 0) & (d < 0.26))[0]:
            if atoms[i].element.symbol == 'O':
                fb.addBond(int(m), int(i), [0.215]); shell0.append(int(i))
    fb.setForceGroup(30)
    system.addForce(fb)
    chi = mm.CustomTorsionForce('0.5*kchi*step(abs(dt)-wchi)*(abs(dt)-wchi)^2; dt = atan2(sin(theta-tchi), cos(theta-tchi))')
    chi.addGlobalParameter('kchi', 1000.0); chi.addGlobalParameter('tchi', math.radians(-124.0))
    chi.addGlobalParameter('wchi', math.radians(35.0))
    for r in top.residues():
        A = {x.name: x.index for x in r.atoms()}
        if r.name not in {'GLY'} | SOLVENT and all(k in A for k in ('N', 'CA', 'C', 'CB')):
            chi.addTorsion(A['N'], A['CA'], A['C'], A['CB'], [])
    chi.setForceGroup(29)
    system.addForce(chi)
    plat = mm.Platform.getPlatformByName(a.platform)
    integ = mm.VerletIntegrator(0.001)
    ctx = mm.Context(system, integ, plat)
    ctx.setPositions(pos)
    ctx.setPeriodicBoxVectors(*top.getPeriodicBoxVectors())
    rec = {'build_dir': os.path.abspath(a.wd), 'platform': a.platform, 'openmm': mm.__version__,
           'n_atoms': system.getNumParticles(), 'mg_shell_start': mg_shell(top, pos)}
    chi_E = lambda: ctx.getState(getEnergy=True, groups={29}).getPotentialEnergy().value_in_unit(KJ)
    ctx.setParameter('kbb', 0.0); ctx.setParameter('kmg', 0.0)
    rec['chirality_restraint_kJmol_start'] = chi_E()
    assert rec['chirality_restraint_kJmol_start'] < 1.0 * max(1, chi.getNumTorsions()) * 0.01 + 5.0, rec['chirality_restraint_kJmol_start']
    ctx.setParameter('kchi', 0.0)
    rec['initial'] = snapshot(ctx, system, top, solute, 'as built (unminimised)')
    ctx.setParameter('kbb', 1000.0); ctx.setParameter('kmg', 50000.0); ctx.setParameter('kchi', 1000.0)
    mm.LocalEnergyMinimizer.minimize(ctx, 100.0, 5000)
    ctx.setParameter('kbb', 0.0); ctx.setParameter('kmg', 0.0); ctx.setParameter('kchi', 0.0)
    rec['after_stage1'] = snapshot(ctx, system, top, solute, 'after restrained minimisation')
    ctx.setParameter('kchi', 1000.0)
    mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 20000)
    rec['chirality_restraint_kJmol_end'] = chi_E()
    ctx.setParameter('kchi', 0.0)
    rec['final'] = snapshot(ctx, system, top, solute, 'after minimisation with only the chirality restraint (tol 10 kJ/mol/nm)')
    st = ctx.getState(getPositions=True)
    P1 = st.getPositions(asNumpy=True)
    X0 = np.array(P0)[bb]; X1 = np.array(P1.value_in_unit(NM))[bb]
    rec['backbone_rmsd_A'] = float(np.sqrt(np.mean(np.sum((X1 - X0) ** 2, 1)))) * 10
    rec['mg_shell_final'] = mg_shell(top, P1)
    if a.chains_a and a.chains_b:
        ca = set(a.chains_a.split(',')); cb = set(a.chains_b.split(','))
        A = [at.index for at in atoms if at.residue.chain.id in ca and at.element.symbol != 'H' and at.residue.name not in SOLVENT]
        B = [at.index for at in atoms if at.residue.chain.id in cb and at.element.symbol != 'H' and at.residue.name not in SOLVENT]
        def contacts(X):
            XA, XB = X[A], X[B]
            n = 0; resA = set(); resB = set()
            for k in range(0, len(XA), 2000):
                d = np.linalg.norm(XA[k:k + 2000, None, :] - XB[None, :, :], axis=2)
                ii, jj = np.where(d < 0.45)
                n += len(ii)
                resA |= {atoms[A[k + i]].residue.index for i in ii}; resB |= {atoms[B[j]].residue.index for j in jj}
            return n, resA, resB
        n0, rA0, rB0 = contacts(np.array(P0)); n1, rA1, rB1 = contacts(np.array(P1.value_in_unit(NM)))
        iface = [at.index for at in atoms if at.residue.index in (rA0 | rB0) and at.name == 'CA']
        d_if = np.sqrt(np.mean(np.sum((np.array(P1.value_in_unit(NM))[iface] - np.array(P0)[iface]) ** 2, 1))) * 10
        rec['interface'] = {'heavy_atom_contacts_4.5A_start': n0, 'heavy_atom_contacts_4.5A_final': n1,
                            'interface_residues_start': len(rA0) + len(rB0), 'interface_residues_final': len(rA1) + len(rB1),
                            'interface_CA_rmsd_A': float(d_if)}
    hand1 = handedness(top, np.array(P1.value_in_unit(NM)))
    flipped = [f'{k[0]}:{k[2]}{k[1]}' for k in hand0 if hand0[k] != hand1.get(k, hand0[k])]
    rec['handedness_changed'] = flipped
    if flipped:
        json.dump(rec, open(os.path.join(a.wd, 'minimize_failed.json'), 'w'), indent=1)
        raise SystemExit(f'handedness changed during minimisation at {flipped}; minimized.cif.gz not written')
    with gzip.open(os.path.join(a.wd, 'minimized.cif.gz'), 'wt') as fh:
        app.PDBxFile.writeFile(top, P1, fh, keepIds=True)
    rec['wall_s'] = round(time.time() - t0, 1)
    json.dump(rec, open(os.path.join(a.wd, 'minimize.json'), 'w'), indent=1)
    for k in ('initial', 'after_stage1', 'final'):
        r = rec[k]
        print(f"{k:13s} E={r['potential_kcalmol']:14.1f} kcal/mol  Fmax={r['max_force_kcalmol_A']:10.2f} kcal/mol/A ({r['max_force_atom']})  "
              f"RMS(all)={r['rms_force_all_kcalmol_A']:.3f}  RMS(solute)={r['rms_force_solute_kcalmol_A']:.3f}")
    print('backbone RMSD', round(rec['backbone_rmsd_A'], 3), 'A | Mg shell final', rec['mg_shell_final'])
    if 'interface' in rec:
        print('interface', rec['interface'])


if __name__ == '__main__':
    main()
