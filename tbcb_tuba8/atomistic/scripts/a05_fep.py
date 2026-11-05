"""Alchemical point-mutation free energy for folding stability, ddG_fold = dG(folded) - dG(unfolded reference).

Hybrid: dual topology from CB outward.  The wild-type side chain (CB and beyond, with H) and the mutant side chain
are both bonded to the same CA; each is its own openmmtools alchemical region (no interactions between them).  The
backbone N,H,CA,HA,C,O is shared.  ff14SB gives CA/HA slightly different charges in different residues; the
difference is added to the mutant CB so that both end states carry integer residue charge (identical in the folded
and the unfolded leg, so it cancels in ddG to first order).  Torsions that involve a side-chain atom are scaled with
that side chain's lambda_torsions, so at either end state the dummy group is held only by its bond and angles to CA
and does not bias phi/psi.

Lambda path (22 states): mutant torsions on; wild-type charges off; sterics swapped (wild type off, mutant on, softcore
alpha 0.5); mutant charges on; wild-type torsions off.  Hamiltonian replica exchange between neighbouring states,
1 ps per iteration, energies of every replica in every state stored for MBAR (a09_fep_analysis.py).

Unfolded reference: Ace-(i-2..i+2)-Nme pentapeptide of the same protein, taken from the minimised build, solvated.

Usage: a05_fep.py <build_dir> <run_dir> <rep> --res 114 --to SER --leg folded|unfolded --ns 5 --budget 2800
"""
import os, sys, json, time, argparse, math
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
from hrex import HREX

BB = {'N', 'H', 'CA', 'HA', 'C', 'O', 'OXT', 'H1', 'H2', 'H3'}
AA1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
       'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
       'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def find_res(top, rid):
    hits = [r for r in top.residues() if r.name not in C.SOLVENT | {'GTP', 'MG'} and int(r.id) == rid]
    assert len(hits) == 1, f'residue {rid}: {len(hits)} matches'
    return hits[0]


def kabsch(X, Y):
    xc, yc = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - xc).T @ (Y - yc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return lambda Z: (R @ (np.asarray(Z) - xc).T).T + yc


def donor_sidechain(top, P, target, newname, exclude_idx, lib=None):
    """Copy a side chain of residue type newname from elsewhere in the same protein (correct chirality, a real
    rotamer, hydrogens included), superposed on the target's N/CA/C; choose the copy with the fewest clashes."""
    env = np.array([a.index for a in top.atoms() if a.element is not None and a.element.symbol != 'H'
                    and a.index not in exclude_idx])
    tgt = {a.name: a.index for a in target.atoms()}
    best = None
    ltop, LP = lib if lib is not None else (top, P)
    for r in ltop.residues():
        if r.name != newname or (lib is None and r.index == target.index):
            continue
        at = {a.name: a.index for a in r.atoms()}
        if not all(k in at for k in ('N', 'CA', 'C')):
            continue
        f = kabsch(LP[[at['N'], at['CA'], at['C']]], P[[tgt['N'], tgt['CA'], tgt['C']]])
        sc = [a for a in r.atoms() if a.name not in BB]
        X = f(LP[[a.index for a in sc]])
        heavy = np.array([a.element.symbol != 'H' for a in sc])
        d = np.linalg.norm(X[heavy][:, None] - P[env][None], axis=2)
        score = float(np.sum(np.clip(0.30 - d, 0, None) ** 2))
        if best is None or score < best[0]:
            bonds = [(b[0].name, b[1].name) for b in ltop.bonds() if b[0].residue is r and b[1].residue is r
                     and (b[0].name not in BB or b[1].name not in BB)]
            best = (score, [(a.name, a.element) for a in sc], X, bonds, int(r.id))
    assert best is not None, f'no donor residue of type {newname}'
    return best


def build_mutant(top, P, target, newname, lib=None):
    """Return mutant topology, positions, map A-index -> B-index, list of new B indices, list of old A side-chain indices."""
    old_sc = [a.index for a in target.atoms() if a.name not in BB]
    score, sc_atoms, X, sc_bonds, donor = donor_sidechain(top, P, target, newname, set(a.index for a in target.atoms()), lib)
    T = app.Topology(); T.setPeriodicBoxVectors(top.getPeriodicBoxVectors())
    amap, new_idx, posB, byname = {}, [], [], {}
    for ch in top.chains():
        c = T.addChain(ch.id)
        for r in ch.residues():
            rr = T.addResidue(newname if r.index == target.index else r.name, c, r.id, r.insertionCode)
            for a in r.atoms():
                if r.index == target.index and a.name not in BB:
                    continue
                b = T.addAtom(a.name, a.element, rr, a.id)
                amap[a.index] = b.index; posB.append(P[a.index])
                if r.index == target.index:
                    byname[a.name] = b
            if r.index == target.index:
                for (nm, el), x in zip(sc_atoms, X):
                    b = T.addAtom(nm, el, rr)
                    new_idx.append(b.index); posB.append(x); byname[nm] = b
    Batoms = list(T.atoms())
    for b1, b2 in top.bonds():
        if b1.index in amap and b2.index in amap:
            T.addBond(Batoms[amap[b1.index]], Batoms[amap[b2.index]])
    for n1, n2 in sc_bonds:
        T.addBond(byname[n1], byname[n2])
    return T, np.array(posB), amap, new_idx, old_sc, dict(donor_residue=donor, donor_clash_score=score)


def hybrid_system(SA, SB, amap, new_idx, old_sc, target_core_A, new_names):
    H = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(SA))
    b2h = {b: a for a, b in amap.items()}
    for b in new_idx:
        b2h[b] = H.addParticle(SB.getParticleMass(b))
    newH = [b2h[b] for b in new_idx]
    newset = set(new_idx)
    for k in range(SB.getNumConstraints()):
        i, j, d = SB.getConstraintParameters(k)
        if i in newset or j in newset:
            H.addConstraint(b2h[i], b2h[j], d)
    get = lambda S, t: [f for f in S.getForces() if type(f).__name__ == t]
    for t in ('HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce', 'NonbondedForce'):
        assert len(get(H, t)) == 1 and len(get(SB, t)) == 1, t
    hb, bb_ = get(H, 'HarmonicBondForce')[0], get(SB, 'HarmonicBondForce')[0]
    for k in range(bb_.getNumBonds()):
        i, j, r0, kb = bb_.getBondParameters(k)
        if i in newset or j in newset:
            hb.addBond(b2h[i], b2h[j], r0, kb)
    ha, ba = get(H, 'HarmonicAngleForce')[0], get(SB, 'HarmonicAngleForce')[0]
    for k in range(ba.getNumAngles()):
        i, j, l, th, ka = ba.getAngleParameters(k)
        if newset & {i, j, l}:
            ha.addAngle(b2h[i], b2h[j], b2h[l], th, ka)
    hp, bp = get(H, 'PeriodicTorsionForce')[0], get(SB, 'PeriodicTorsionForce')[0]
    for k in range(bp.getNumTorsions()):
        i, j, l, m, per, ph, kt = bp.getTorsionParameters(k)
        if newset & {i, j, l, m}:
            hp.addTorsion(b2h[i], b2h[j], b2h[l], b2h[m], per, ph, kt)
    hn, bn = get(H, 'NonbondedForce')[0], get(SB, 'NonbondedForce')[0]
    qA = lambda i: hn.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
    qB = lambda i: bn.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
    # every retained atom must have identical LJ, and identical charge outside the mutated residue's backbone
    core = set(target_core_A)
    for a, b in amap.items():
        pa, pb = hn.getParticleParameters(a), bn.getParticleParameters(b)
        assert abs(pa[1].value_in_unit(unit.nanometer) - pb[1].value_in_unit(unit.nanometer)) < 1e-6
        assert abs(pa[2].value_in_unit(unit.kilojoule_per_mole) - pb[2].value_in_unit(unit.kilojoule_per_mole)) < 1e-6
        if a not in core:
            assert abs(qA(a) - qB(b)) < 1e-6, (a, qA(a), qB(b))
    dq = sum(qB(amap[a]) - qA(a) for a in target_core_A)
    qH = {}
    for b, nm in zip(new_idx, new_names):
        q, s, e = bn.getParticleParameters(b)
        qv = q.value_in_unit(unit.elementary_charge) + (dq if nm == 'CB' else 0.0)
        hn.addParticle(qv, s, e); qH[b2h[b]] = qv
    for i in range(H.getNumParticles()):
        if i not in qH:
            qH[i] = qA(i)
    scale14 = []
    for k in range(bn.getNumExceptions()):
        i, j, qq, s, e = bn.getExceptionParameters(k)
        if not (i in newset or j in newset):
            continue
        qq = qq.value_in_unit(unit.elementary_charge ** 2)
        if qq != 0.0:
            scale14.append(qq / (qB(i) * qB(j)))
            qq = scale14[-1] * qH[b2h[i]] * qH[b2h[j]]
        hn.addException(b2h[i], b2h[j], qq, s, e)
    for i in old_sc:
        for j in newH:
            hn.addException(i, j, 0.0, 0.1, 0.0)
    s14 = float(np.median(scale14)) if scale14 else None
    assert s14 is None or abs(s14 - 1 / 1.2) < 1e-3, s14
    cbB = [b for b, nm in zip(new_idx, new_names) if nm == 'CB']
    return H, newH, dict(core_charge_shift_added_to_CB=dq, scale14=s14, _cb_B=cbB[0] if cbB else -1)


def pentapeptide(top, P, rid):
    """Ace-(rid-2..rid+2)-Nme from the folded structure; caps are built from residues rid-3 and rid+3 by turning their
    CA into the methyl carbon and placing methyl H along the old CA bonds."""
    res = {int(r.id): r for r in top.residues() if r.name not in C.SOLVENT | {'GTP', 'MG'}}
    T = app.Topology(); ch = T.addChain('P'); pos = []
    def at(r, n):
        return [a for a in r.atoms() if a.name == n][0]
    def methyl_H(cpos, partners):
        out = []
        for p in partners:
            v = p - cpos; out.append(cpos + 0.109 * v / np.linalg.norm(v))
        return out
    H = app.element.hydrogen; Cel = app.element.carbon; Oel = app.element.oxygen; Nel = app.element.nitrogen
    rA = res[rid - 3]; ace = T.addResidue('ACE', ch, str(rid - 3))
    ca = P[at(rA, 'CA').index]
    partners = [P[at(rA, 'N').index]] + [P[a.index] for a in rA.atoms() if a.name in ('HA', 'CB', 'HA2', 'HA3')]
    hs = methyl_H(ca, partners[:3])
    aC = T.addAtom('C', Cel, ace); pos.append(P[at(rA, 'C').index])
    aO = T.addAtom('O', Oel, ace); pos.append(P[at(rA, 'O').index])
    aM = T.addAtom('CH3', Cel, ace); pos.append(ca)
    hh = []
    for k, x in enumerate(hs):
        hh.append(T.addAtom(f'HH3{k + 1}', H, ace)); pos.append(x)
    T.addBond(aC, aO); T.addBond(aC, aM); [T.addBond(aM, h) for h in hh]
    prevC = aC; keep = {}
    for k in range(rid - 2, rid + 3):
        r = res[k]; rr = T.addResidue(r.name, ch, r.id)
        for a in r.atoms():
            keep[a.index] = T.addAtom(a.name, a.element, rr); pos.append(P[a.index])
        T.addBond(prevC, keep[at(r, 'N').index]); prevC = keep[at(r, 'C').index]
    for b1, b2 in top.bonds():
        if b1.index in keep and b2.index in keep:
            T.addBond(keep[b1.index], keep[b2.index])
    rN = res[rid + 3]; nme = T.addResidue('NME', ch, str(rid + 3))
    nN = T.addAtom('N', Nel, nme); pos.append(P[at(rN, 'N').index])
    nH = T.addAtom('H', H, nme); pos.append(P[at(rN, 'H').index])
    ca = P[at(rN, 'CA').index]
    partners = [P[at(rN, 'C').index]] + [P[a.index] for a in rN.atoms() if a.name in ('HA', 'CB', 'HA2', 'HA3')]
    nM = T.addAtom('CH3', Cel, nme); pos.append(ca)
    T.addBond(prevC, nN); T.addBond(nN, nH); T.addBond(nN, nM)
    for k, x in enumerate(methyl_H(ca, partners[:3])):
        h = T.addAtom(f'HH3{k + 1}', H, nme); pos.append(x); T.addBond(nM, h)
    return T, np.array(pos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('build'); ap.add_argument('run'); ap.add_argument('rep', type=int)
    ap.add_argument('--res', type=int, required=True); ap.add_argument('--to', required=True)
    ap.add_argument('--leg', choices=['folded', 'unfolded'], required=True)
    ap.add_argument('--ns', type=float, default=5.0); ap.add_argument('--budget', type=float, default=2800.0)
    ap.add_argument('--steps', type=int, default=250); ap.add_argument('--dt', type=float, default=0.004)
    ap.add_argument('--platform', default=None); ap.add_argument('--eq_scale', type=float, default=1.0)
    ap.add_argument('--test', action='store_true', help='endpoint consistency test only')
    ap.add_argument('--donor_lib', default=None, help='structure to copy the new side chain from (default: the build)')
    ap.add_argument('--restrain', default=None,
                    help='flat-bottom pose restraint "A:187,203,204,223|B:451": CA centroid of the first set vs N/CA/C '
                         'centroid of the second, inactive below start distance + 0.3 nm (same in every state)')
    a = ap.parse_args()
    W = os.path.join(a.run, f'r{a.rep}'); os.makedirs(W, exist_ok=True)
    LOCAL = os.path.join(os.environ.get('L_SCRATCH_JOB') or os.environ.get('TMPDIR') or '/tmp', f'fep_{os.path.basename(a.run)}_r{a.rep}')
    os.makedirs(LOCAL, exist_ok=True)
    log = lambda m: print(m, flush=True)
    ff = C.forcefield()
    seed = 777 + 31 * a.rep + 7 * (a.leg == 'unfolded')
    EQ = os.path.join(W, 'equil_state.npz'); TOPF = os.path.join(W, 'wt_system.cif.gz')
    if not os.path.exists(EQ) or a.test:
        top0, pos0 = C.load_cif(os.path.join(a.build, 'minimized.cif.gz'))
        if a.leg == 'unfolded':
            T, X = pentapeptide(top0, np.array(pos0.value_in_unit(C.NM)), a.res)
            m = app.Modeller(T, X * C.NM)
            # cubic box, never narrower than 3.6 nm, so NPT fluctuations cannot bring it near twice the 1.0 nm cutoff
            ext = float(np.max(X.max(0) - X.min(0)))
            L = max(ext + 2.4, 3.6)
            m.addSolvent(ff, model='tip3p', boxSize=mm.Vec3(L, L, L) * C.NM, ionicStrength=0.15 * unit.molar,
                         neutralize=True)
            S = C.make_system(m.topology, hmr=False, constraints=None)
            ctx = mm.Context(S, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference' if S.getNumParticles() < 200 else 'CPU'))
            ctx.setPositions(m.positions); mm.LocalEnergyMinimizer.minimize(ctx, 10.0, 5000)
            top, pos = m.topology, ctx.getState(getPositions=True).getPositions()
        else:
            top, pos = top0, pos0
        C.write_cif(TOPF, top, pos)
        if a.test:
            P = np.array(pos.value_in_unit(C.NM)); V = np.zeros_like(P); box = top.getPeriodicBoxVectors()
            C.save_state_npz(os.path.join(W, 'test_state.npz'), P, V, box)
        else:
            P, V, box, rec = C.equilibrate(top, pos, T=310.0, seed=seed, platform=a.platform, log=log, npt=True, scale=a.eq_scale)
            json.dump(rec, open(os.path.join(W, 'equil.json'), 'w'), indent=1)
            C.save_state_npz(EQ, P, V, box)
    top, _ = C.load_cif(TOPF)
    P, V, box, _ = C.load_state_npz(os.path.join(W, 'test_state.npz') if a.test else EQ)
    top.setPeriodicBoxVectors(box)
    target = find_res(top, a.res)
    oldname = target.name
    assert oldname != 'GLY' and a.to != 'GLY' and 'PRO' not in (oldname, a.to)
    ltop, lpos = C.load_cif(a.donor_lib or os.path.join(a.build, 'minimized.cif.gz'))
    TB, PB, amap, new_idx, old_sc, dinfo = build_mutant(top, P, target, a.to, lib=(ltop, np.array(lpos.value_in_unit(C.NM))))
    hmr = a.dt > 0.0025
    SA = C.make_system(top, hmr=hmr, ff=ff); SB = C.make_system(TB, hmr=hmr, ff=ff)
    for S in (SA, SB):
        [f for f in S.getForces() if isinstance(f, mm.NonbondedForce)][0].setUseDispersionCorrection(False)
    core_A = [x.index for x in target.atoms() if x.name in BB]
    Batoms = list(TB.atoms())
    new_names = [Batoms[b].name for b in new_idx]
    H, newH, hinfo = hybrid_system(SA, SB, amap, new_idx, old_sc, core_A, new_names)
    from openmmtools import alchemy
    # Aromatic side chains: switching ALL torsions of a dummy ring makes the ring floppy at the dummy end, and turning
    # its ~3.6 kcal/mol ring torsions back on gives near-zero exchange between the first/last lambda pairs (seen in
    # the first Y<->F runs: acceptance 0.00-0.04).  For aromatic mutations only the torsions that cross from the
    # side chain into the rest of the system are scaled; torsions internal to the side chain stay on, which keeps
    # the dummy ring planar.  Internal terms of a non-interacting group factorise and cancel between legs.
    AROM = {'PHE', 'TYR', 'TRP', 'HIS', 'HIE', 'HID', 'HIP'}
    tor_old, tor_new = True, True
    if oldname in AROM or a.to in AROM:
        ptH = [f for f in H.getForces() if isinstance(f, mm.PeriodicTorsionForce)][0]
        oset, nset = set(old_sc), set(newH)
        tor_old, tor_new = [], []
        for k in range(ptH.getNumTorsions()):
            i_, j_, l_, m_, _, _, _ = ptH.getTorsionParameters(k)
            at4 = {i_, j_, l_, m_}
            if at4 & oset and at4 - oset:
                tor_old.append(k)
            if at4 & nset and at4 - nset:
                tor_new.append(k)
    reg_old = alchemy.AlchemicalRegion(alchemical_atoms=old_sc, alchemical_torsions=tor_old, annihilate_electrostatics=True,
                                       annihilate_sterics=False, softcore_alpha=0.5, name='old')
    reg_new = alchemy.AlchemicalRegion(alchemical_atoms=newH, alchemical_torsions=tor_new, annihilate_electrostatics=True,
                                       annihilate_sterics=False, softcore_alpha=0.5, name='new')
    hinfo['torsion_scaling'] = 'all side-chain torsions' if tor_old is True else f'boundary only ({len(tor_old)} old, {len(tor_new)} new)'
    fac = alchemy.AbsoluteAlchemicalFactory(alchemical_pme_treatment='exact', disable_alchemical_dispersion_correction=True)
    AS = fac.create_alchemical_system(H, [reg_old, reg_new])
    cvf = None; rinfo = None
    if a.restrain:
        def sel(spec, names):
            cid, ids = spec.split(':'); ids = {int(x) for x in ids.split(',')}
            return [x.index for x in top.atoms() if x.residue.chain.id == cid and int(x.residue.id) in ids and x.name in names]
        sA, sB = a.restrain.split('|')
        g1 = sel(sA, ('CA',)); g2 = sel(sB, ('N', 'CA', 'C'))
        assert g1 and g2, (g1, g2)
        Pnm = np.asarray(P)
        d_init = float(np.linalg.norm(Pnm[g1].mean(0) - Pnm[g2].mean(0)))
        fb = mm.CustomCentroidBondForce(2, '0.5*kpose*step(distance(g1,g2)-d0)*(distance(g1,g2)-d0)^2')
        fb.addGlobalParameter('kpose', 1000.0); fb.addPerBondParameter('d0')
        fb.addGroup(g1); fb.addGroup(g2); fb.addBond([0, 1], [d_init + 0.3])
        fb.setUsesPeriodicBoundaryConditions(True)
        AS.addForce(fb)
        rinfo = dict(restraint_groups=[len(g1), len(g2)], d_start_nm=d_init, d0_nm=d_init + 0.3, k=1000.0)
        def cvf(X):
            return [float(np.linalg.norm(X[g1].mean(0) - X[g2].mean(0)))]
    XH = np.concatenate([P, PB[new_idx]], axis=0)
    probe = mm.Context(AS, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU' if a.test else (a.platform or 'CUDA')))
    gp = sorted(probe.getParameters().keys()) if hasattr(probe.getParameters(), 'keys') else []
    lam_params = [p for p in gp if p.startswith('lambda_')]
    log(f'mutation {AA1[oldname]}{a.res}{AA1[a.to]} leg {a.leg}: {len(old_sc)} old / {len(newH)} new atoms; '
        f'parameters {lam_params}; {json.dumps(dinfo)}; {json.dumps(hinfo)}')
    del probe
    base = {p: 1.0 for p in lam_params}
    def st(eo, so, to, en, sn, tn):
        s = dict(base)
        s.update({'lambda_electrostatics_old': eo, 'lambda_sterics_old': so, 'lambda_torsions_old': to,
                  'lambda_electrostatics_new': en, 'lambda_sterics_new': sn, 'lambda_torsions_new': tn})
        return s
    states = [st(1, 1, 1, 0, 0, 0), st(1, 1, 1, 0, 0, .5), st(1, 1, 1, 0, 0, 1)]
    states += [st(e, 1, 1, 0, 0, 1) for e in (.75, .5, .25, 0)]
    arom_swap = oldname in AROM or a.to in AROM
    xs = (.05, .1, .15, .2, .25, .3, .35, .45, .55, .7, .85, .95) if arom_swap else (.05, .15, .3, .45, .55, .7, .85, .95)
    states += [st(0, 1 - x, 1, 0, x, 1) for x in xs]
    states += [st(0, 0, 1, 0, 1, 1)]
    states += [st(0, 0, 1, e, 1, 1) for e in (.25, .5, .75, 1)]
    states += [st(0, 0, .5, 1, 1, 1), st(0, 0, 0, 1, 1, 1)]
    for s in states:
        assert set(s) == set(lam_params), (sorted(s), lam_params)
    if a.test:
        endpoint_test(SA, SB, AS, top, P, PB, amap, new_idx, old_sc, newH, states, box, log, W, target, hinfo)
        return
    meta = dict(mutation=f'{AA1[oldname]}{a.res}{AA1[a.to]}', leg=a.leg, build=a.build, n_states=len(states), states=states,
                n_atoms=AS.getNumParticles(), old_sc=len(old_sc), new_sc=len(newH), dt_ps=a.dt, steps_per_iter=a.steps,
                pose_restraint=rinfo, **dinfo, **hinfo)
    json.dump(meta, open(os.path.join(W, 'fep_meta.json'), 'w'), indent=1)
    n_iter = int(round(a.ns * 1000 / (a.steps * a.dt)))
    h = HREX(AS, XH, box, states, temperature_K=310.0, timestep_ps=a.dt, steps_per_iter=a.steps, platform=a.platform,
             shared_dir=W, local_dir=LOCAL, tag='fep', seed=seed, cv_fn=cvf)
    # walk the ladder from state 0 so every state starts from a configuration relaxed in its own Hamiltonian
    h.resume_or_init(sequential=True, warm_steps=2500)
    done = h.run(n_iter, a.budget, ck_every_min=10.0, log_every=100)
    json.dump({'iteration': h.iteration, 'target': n_iter, 'done': bool(done), 'acceptance': h.acceptance(),
               'platform': h.platform, 'time': time.ctime()}, open(os.path.join(W, 'fep_status.json'), 'w'), indent=1)
    log(f'iteration {h.iteration}/{n_iter} done={done} acceptance min {min(h.acceptance()):.2f}')


def endpoint_test(SA, SB, AS, top, P, PB, amap, new_idx, old_sc, newH, states, box, log, W, target, hinfo):
    """Endpoints must reproduce the physical wild-type and mutant energies up to a constant that depends only on the
    dummy group and the mutated residue (its own bonded and intra-group terms).  Perturb every atom outside the
    mutated residue and check that hybrid - physical is unchanged.  The mutant reference carries the same
    CA/HA -> CB charge transfer as the hybrid (otherwise the check would be inexact by construction)."""
    plat = mm.Platform.getPlatformByName('CPU')
    SBm = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(SB))
    nbB = [f for f in SBm.getForces() if isinstance(f, mm.NonbondedForce)][0]
    nbA = [f for f in SA.getForces() if isinstance(f, mm.NonbondedForce)][0]
    core = [x.index for x in target.atoms() if x.name in BB]
    Batoms = list(SBm.getNumParticles() * [None])
    qnew = {}
    for i in core:
        b = amap[i]; q, s_, e = nbB.getParticleParameters(b); qa = nbA.getParticleParameters(i)[0]
        nbB.setParticleParameters(b, qa, s_, e); qnew[b] = qa.value_in_unit(unit.elementary_charge)
    names = {}
    for ch in []:
        pass
    cb = [b for b in new_idx if b in hinfo.get('_cb', [])]
    return_cb = None
    for b in new_idx:
        q, s_, e = nbB.getParticleParameters(b)
        if abs(hinfo['core_charge_shift_added_to_CB']) > 0 and b == hinfo['_cb_B']:
            q = q.value_in_unit(unit.elementary_charge) + hinfo['core_charge_shift_added_to_CB']
            nbB.setParticleParameters(b, q, s_, e); qnew[b] = q
    qq = lambda i: nbB.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
    for k in range(nbB.getNumExceptions()):
        i, j, c, s_, e = nbB.getExceptionParameters(k)
        if (i in qnew or j in qnew) and c.value_in_unit(unit.elementary_charge ** 2) != 0:
            nbB.setExceptionParameters(k, i, j, hinfo['scale14'] * qq(i) * qq(j), s_, e)
    def E(S, X, params=None):
        c = mm.Context(S, mm.VerletIntegrator(0.001), plat)
        c.setPeriodicBoxVectors(*box); c.setPositions(X)
        if params:
            for k_, v in params.items():
                c.setParameter(k_, v)
        return c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(C.KJ)
    rng = np.random.default_rng(0)
    nA = SA.getNumParticles()
    fixed = set(x.index for x in target.atoms())
    mask = np.array([i not in fixed for i in range(nA)])
    XB_new = PB[new_idx]
    out = {}
    for trial in range(3):
        XA = P.copy()
        if trial:
            XA[mask] += rng.normal(0, 0.002, size=(mask.sum(), 3))
        XH = np.concatenate([XA, XB_new])
        XB = np.zeros((SB.getNumParticles(), 3))
        for i, b in amap.items():
            XB[b] = XA[i]
        XB[new_idx] = XB_new
        eA = E(SA, XA); eB = E(SBm, XB)
        h0 = E(AS, XH, states[0]); hK = E(AS, XH, states[-1])
        out[trial] = dict(hybrid0_minus_WT=h0 - eA, hybridK_minus_MUT=hK - eB, E_WT=eA, E_MUT=eB)
        log(f'trial {trial}: hybrid(state 0) - WT = {h0 - eA:.5f}   hybrid(state K) - MUT = {hK - eB:.5f} kJ/mol')
    drift = {k: max(abs(out[t][k] - out[0][k]) for t in (1, 2)) for k in ('hybrid0_minus_WT', 'hybridK_minus_MUT')}
    log(f'environment-perturbation drift, must be ~0: {drift}')
    json.dump({'trials': out, 'max_drift_kJmol': drift}, open(os.path.join(W, 'endpoint_test.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
