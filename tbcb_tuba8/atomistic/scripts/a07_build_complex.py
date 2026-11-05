"""Build solvated complexes that contain alpha-tubulin, keeping GTP.Mg2+ at its N-site.

Sources:
  af3   : the user's AlphaFold3 models (no ligands).  GTP + Mg2+ + its four inner-shell waters are grafted from the
          1.8 A crystal 4I4T chain A, superposed on the alpha N-site residues.  Clashes of the graft with the partner
          are recorded (the AF3 TBCB pose sits on the beta-binding face, near the N-site).
  exp   : an experimental cofactor complex (9M1L / 9M1M), chains given explicitly; alpha can be re-threaded to another
          isoform exactly as in a01_build.py; GTP/Mg of every tubulin chain kept, Mg shells completed from 4I4T.
The solute is rotated to its principal axes and solvated in a rectangular box (1.2 nm padding, 0.15 M NaCl).

Usage:
  a07_build_complex.py <sid> af3 <model.cif> --alpha B --partners A
  a07_build_complex.py <sid> exp <9M1L.cif.gz> --alpha a --partners E [--thread Tuba8] [--keep-beta b]
Writes build/<sid>/{solvated.cif.gz, build.json, start_noH.pdb, vacuum_forces.json}
"""
import os, sys, io, json, gzip, argparse
import numpy as np
import gemmi
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a01_build as A1

ROOT = A1.ROOT; FF = A1.FF


def kabsch(X, Y):
    xc, yc = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - xc).T @ (Y - yc)); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return lambda Z: (R @ (np.asarray(Z) - xc).T).T + yc


def graft_nsite(alpha_chain, ref_chain_name='A', ref=os.path.join(ROOT, 'inputs', '4I4T.cif.gz')):
    """GTP, Mg and 4 Mg-bound waters from 4I4T placed on alpha_chain (gemmi chain).  Returns PDB HETATM lines and info."""
    rs = gemmi.read_structure(ref); rs.setup_entities(); rch = rs[0][ref_chain_name]
    rg = [r for r in rch if r.name == 'GTP'][0]; rmg = [r for r in rch if r.name == 'MG'][0]
    rgp = np.array([a.pos.tolist() for a in rg if a.element.name != 'H'])
    site = sorted({r.seqid.num for r in rch.get_polymer() for a in r
                   if np.min(np.linalg.norm(rgp - np.array(a.pos.tolist()), axis=1)) < 6.0})
    tca = {r.seqid.num: np.array(r['CA'][0].pos.tolist()) for r in alpha_chain if r.find_atom('CA', '*')}
    rca = {r.seqid.num: np.array(r['CA'][0].pos.tolist()) for r in rch.get_polymer() if r.find_atom('CA', '*')}
    common = [p for p in site if p in tca and p in rca]
    f = kabsch(np.array([rca[p] for p in common]), np.array([tca[p] for p in common]))
    rmsd = float(np.sqrt(np.mean(np.sum((f(np.array([rca[p] for p in common])) - np.array([tca[p] for p in common])) ** 2, 1))))
    ns = gemmi.NeighborSearch(rs[0], rs.cell, 4).populate()
    shell = [f(mk.to_cra(rs[0]).atom.pos.tolist()) for mk in ns.find_atoms(rmg[0].pos, '\0', radius=2.6)
             if mk.to_cra(rs[0]).residue.name == 'HOH']
    assert len(shell) == 4
    mgpos = f(rmg[0].pos.tolist())
    atoms = [(a.name, a.element.name, f(a.pos.tolist())) for a in rg if a.element.name != 'H']
    return atoms, mgpos, shell, dict(site_ca_rmsd_A=round(rmsd, 3), n_site=len(common))


def het_lines(atoms, mgpos, shell, chain_gtp='X', chain_mg='Y', chain_w='Z', start=1, resnum=901):
    lines = []; k = start
    for nm, el, x in atoms:
        lines.append(f"HETATM{k:5d} {nm:<4s} GTP {chain_gtp}{resnum:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 20.00          {el:>2s}"); k += 1
    lines.append(f"HETATM{k:5d} MG    MG {chain_mg}{resnum:4d}    {mgpos[0]:8.3f}{mgpos[1]:8.3f}{mgpos[2]:8.3f}  1.00 20.00          MG"); k += 1
    for i, o in enumerate(shell):
        for nm, x in zip(('O', 'H1', 'H2'), A1.water_atoms(o, mgpos)):
            lines.append(f"HETATM{k:5d} {nm:<4s} HOH {chain_w}{resnum + i:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 20.00           {nm[0]}"); k += 1
    return lines, k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sid'); ap.add_argument('kind', choices=['af3', 'exp']); ap.add_argument('src')
    ap.add_argument('--alpha', required=True); ap.add_argument('--partners', required=True)
    ap.add_argument('--thread', default=None); ap.add_argument('--keep-beta', default=None)
    ap.add_argument('--pad', type=float, default=1.2); ap.add_argument('--out', default=os.path.join(ROOT, 'build'))
    ap.add_argument('--vacuum-only', action='store_true')
    a = ap.parse_args()
    W = os.path.join(a.out, a.sid); os.makedirs(W, exist_ok=True)
    st = gemmi.read_structure(a.src); st.setup_entities(); st.remove_hydrogens()
    keep = list(dict.fromkeys([a.alpha] + a.partners.split(',') + ([a.keep_beta] if a.keep_beta else [])))
    info = dict(sid=a.sid, kind=a.kind, src=os.path.basename(a.src), alpha=a.alpha, partners=a.partners, beta=a.keep_beta)
    aln = A1.load_alignment()
    if a.kind == 'af3':
        out = gemmi.Structure(); out.add_model(gemmi.Model('1'))
        for cid in keep:
            ch = st[0][cid]; nc = gemmi.Chain(cid)
            for r in ch.get_polymer():
                nc.add_residue(r)
            out[0].add_chain(nc)
        out.setup_entities()
        fixer = PDBFixer(pdbfile=io.StringIO(out.make_pdb_string()))
        fixer.missingResidues = {}
    else:
        # experimental: PDBFixer straight from the mmCIF so SEQRES is known; build internal gaps only (never termini,
        # so the unresolved alpha C-terminal tail stays absent as in the monomers); thread alpha by applyMutations
        fixer = PDBFixer(filename=a.src)
        # PDBFixer names chains by label_asym_id, gemmi by auth; match them by residue count and first residue
        lab = {}
        for cid in keep:
            g = st[0][cid].get_polymer()
            sig = (len(g), g[0].name, str(g[0].seqid.num))
            hit = [c.id for c in fixer.topology.chains()
                   if (sum(1 for _ in c.residues()), next(c.residues()).name, next(c.residues()).id) == sig]
            assert len(hit) == 1, (cid, sig, hit)
            lab[hit[0]] = cid
        ids = [c.id for c in fixer.topology.chains()]
        fixer.removeChains(chainIds=[c for c in ids if c not in lab])
        fixer.findMissingResidues()   # must run before the chains are renamed: PDBFixer matches SEQRES by chain id
        for c in fixer.topology.chains():
            c.id = lab[c.id]
        info['label_to_auth'] = lab
        chains = list(fixer.topology.chains())
        internal = {}
        for (ci, ri), resl in fixer.missingResidues.items():
            n = len(list(chains[ci].residues()))
            if 0 < ri < n:
                internal[(ci, ri)] = resl
        fixer.missingResidues = internal
        info['internal_gaps_built'] = {f'{chains[ci].id}:{ri}': len(v) for (ci, ri), v in internal.items()}
        fixer.removeHeterogens(keepWater=False)
        if a.thread:
            muts = []
            for r in [c for c in fixer.topology.chains() if c.id == a.alpha][0].residues():
                want = aln[a.thread].get(int(r.id))
                if want and want != '-' and A1.AA3[want] != r.name:
                    muts.append(f'{r.name}-{r.id}-{A1.AA3[want]}')
            info['threaded_mutations'] = muts
            if muts:
                fixer.applyMutations(muts, a.alpha)
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    app.Topology.loadBondDefinitions(os.path.join(FF, 'gtp_bonds.xml'))
    het = []; k = 1; grafts = {}
    tub = [a.alpha]
    for t_i, cid in enumerate(tub):
        g_atoms, mgpos, shell, gi = graft_nsite(st[0][cid], ref_chain_name='A')
        lines, k = het_lines(g_atoms, mgpos, shell, chain_gtp='XW'[t_i], chain_mg='YV'[t_i], chain_w='ZU'[t_i], start=k,
                             resnum=901 + 10 * t_i)
        het += lines; grafts[cid] = gi
        partner_heavy = np.array([at.pos.tolist() for c in a.partners.split(',') for r in st[0][c].get_polymer() for at in r])
        own_heavy = np.array([at.pos.tolist() for r in st[0][cid].get_polymer() for at in r])
        G = np.array([x for _, _, x in g_atoms] + [mgpos] + shell)
        gi['min_dist_graft_to_partner_A'] = round(float(np.min(np.linalg.norm(G[:, None] - partner_heavy[None], axis=2))), 2)
        gi['min_dist_graft_to_own_chain_A'] = round(float(np.min(np.linalg.norm(G[:, None] - own_heavy[None], axis=2))), 2)
        gi['n_partner_heavy_within_3A_of_graft'] = int(np.sum(np.min(np.linalg.norm(partner_heavy[:, None] - G[None], axis=2), axis=1) < 3.0))
    if a.keep_beta:
        # beta keeps its own GTP + Mg from the experimental structure; the Mg shell is completed with waters on the
        # octahedral sites left free by the phosphate oxygens (a01_build.octahedral_waters)
        bch = st[0][a.keep_beta]
        blig = [r for r in bch.get_ligands() if r.name in ('GTP', 'GDP', 'MG')]
        bg = [r for r in blig if r.name in ('GTP', 'GDP')][0]; bmg = [r for r in blig if r.name == 'MG']
        gat = [(x.name, x.element.name, np.array(x.pos.tolist())) for x in bg if x.element.name != 'H']
        if bmg:
            mgp = np.array(bmg[0][0].pos.tolist())
            near = sorted([(np.linalg.norm(x - mgp), n, x) for n, el, x in gat if el == 'O'], key=lambda t: t[0])[:2]
            protO = np.array([np.array(x.pos.tolist()) for r in bch.get_polymer() for x in r if x.element.name == 'O'])
            ws = A1.octahedral_waters(mgp, [near[0][2], near[1][2]], protO)
            lines, k = het_lines(gat, mgp, ws, chain_gtp='Q', chain_mg='R', chain_w='S', start=k, resnum=951)
        else:
            lines = []
            for nm, el, x in gat:
                lines.append(f"HETATM{k:5d} {nm:<4s} {bg.name} Q 951    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 20.00          {el:>2s}"); k += 1
        het += lines
        grafts[a.keep_beta] = dict(nucleotide=bg.name, mg=bool(bmg), shell_waters_added=len(ws) if bmg else 0)
    lig = app.PDBFile(io.StringIO('\n'.join(het) + '\nEND\n'))
    m = app.Modeller(fixer.topology, fixer.positions)
    m.add(lig.topology, lig.positions)
    with open(os.path.join(W, 'start_noH.pdb'), 'w') as fh:
        app.PDBFile.writeFile(m.topology, m.positions, fh)
    ff = A1.forcefield()
    app.Modeller.loadHydrogenDefinitions(os.path.join(FF, 'gtp_hydrogens.xml'))
    m.addHydrogens(ff, pH=7.0)
    # internal gaps rebuilt by PDBFixer start with overlapping atoms; relax only those residues (+-2) in vacuum with
    # everything else held fixed, before the force record and before solvation
    built = []
    if a.kind == 'exp' and info.get('internal_gaps_built'):
        for key, nres in info['internal_gaps_built'].items():
            cid, ri = key.split(':'); ri = int(ri)
            ch = [c for c in m.topology.chains() if c.id == cid][0]
            rl = list(ch.residues())
            built += [r.index for r in rl[max(0, ri - 2): ri + nres + 2]]
    if built:
        sb = ff.createSystem(m.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        bset = set(built)
        for at in m.topology.atoms():
            if at.residue.index not in bset:
                sb.setParticleMass(at.index, 0.0)
        cb = mm.Context(sb, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
        cb.setPositions(m.positions)
        e0 = cb.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        mm.LocalEnergyMinimizer.minimize(cb, 10.0, 5000)
        e1 = cb.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        m.positions = cb.getState(getPositions=True).getPositions()
        info['built_loop_relaxation'] = dict(n_residues=len(built), E_before_kcal=e0, E_after_kcal=e1)
        del cb
    # vacuum force record of the as-built complex (protein + GTP/Mg, no solvent, no cutoff): the analogue of the
    # "AF3 complex as deposited" number, before any force-field relaxation
    sv = ff.createSystem(m.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    cv = mm.Context(sv, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('CPU'))
    cv.setPositions(m.positions)
    s0 = cv.getState(getEnergy=True, getForces=True)
    Fv = np.linalg.norm(s0.getForces(asNumpy=True).value_in_unit(unit.kilocalorie_per_mole / unit.angstrom), axis=1)
    atoms_ = list(m.topology.atoms()); imax = int(np.argmax(Fv))
    vac = dict(potential_kcalmol=s0.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole),
               rms_force_kcalmol_A=float(np.sqrt(np.mean(Fv ** 2))), max_force_kcalmol_A=float(Fv[imax]),
               max_force_atom=f'{atoms_[imax].residue.chain.id}:{atoms_[imax].residue.name}{atoms_[imax].residue.id}:{atoms_[imax].name}',
               n_atoms_force_gt_100=int(np.sum(Fv > 100)))
    vac['grafts'] = grafts
    json.dump(vac, open(os.path.join(W, 'vacuum_forces.json'), 'w'), indent=1)
    del cv
    if a.vacuum_only:
        print(json.dumps({'sid': a.sid, 'vacuum_as_built': vac})); return
    # principal axes, then a rectangular box
    X = np.array(m.positions.value_in_unit(unit.nanometer))
    c = X.mean(0); w, v = np.linalg.eigh((X - c).T @ (X - c))
    Y = (X - c) @ v
    ext = Y.max(0) - Y.min(0)
    Y = Y - Y.min(0) + a.pad
    box = ext + 2 * a.pad
    m2 = app.Modeller(m.topology, [mm.Vec3(*y) for y in Y] * unit.nanometer)
    m2.addSolvent(ff, model='tip3p', boxSize=mm.Vec3(*box) * unit.nanometer, ionicStrength=0.15 * unit.molar, neutralize=True)
    s = ff.createSystem(m2.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer, constraints=app.HBonds)
    nbf = [f for f in s.getForces() if isinstance(f, mm.NonbondedForce)][0]
    q = sum(nbf.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(s.getNumParticles()))
    with gzip.open(os.path.join(W, 'solvated.cif.gz'), 'wt') as fh:
        app.PDBxFile.writeFile(m2.topology, m2.positions, fh, keepIds=True)
    chains = {ch.id: sum(1 for _ in ch.residues()) for ch in m2.topology.chains() if ch.id in keep}
    info.update(n_atoms=s.getNumParticles(), net_charge=round(q, 6), box_nm=[round(float(b), 3) for b in box],
                chain_residues=chains, grafts=grafts, vacuum_as_built=vac, thread=a.thread, openmm=mm.__version__)
    json.dump(info, open(os.path.join(W, 'build.json'), 'w'), indent=1)
    print(json.dumps({k: info[k] for k in ('sid', 'n_atoms', 'net_charge', 'box_nm', 'grafts', 'vacuum_as_built')}))


if __name__ == '__main__':
    main()
