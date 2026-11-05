"""Build solvated alpha-tubulin systems that keep GTP.Mg2+ at the N-site.

Template: 9M1M chain a (human = mouse TUBA1B, residues 1-437 all resolved, GTP + Mg bound), i.e. alpha-tubulin in a
cofactor-bound state.  Every isoform/mutant is threaded onto the same backbone so that differences between systems
come from sequence, not from the starting model.  Residues that differ from the template keep N/CA/C/O(/CB) and get
their side chain rebuilt by PDBFixer; the unresolved C-terminal tail (438-end) is built by PDBFixer.

Mg2+: 9M1M has no modelled waters and its Mg sits 2.8-3.0 A from O1B, O3G and Glu71 OE2.  The octahedral shell is
completed with explicit waters on the free sites (a site is left empty if a protein oxygen already occupies it), and
the build records the shell.  Coordination is checked again after minimisation.

Usage:  a01_build.py <system_id> <spec> [--pad 1.0]
  spec examples:  Tuba1a     Tuba8:S114I     Tuba8:loop=Tuba1a     Tuba1a:loop=Tuba8,I114S
Outputs (build/<system_id>/): solvated.cif.gz, build.json, start_noH.pdb
"""
import os, sys, json, gzip, csv, io, argparse, math
import numpy as np
import gemmi
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
FF = os.path.join(ROOT, 'ff')
ALN = os.path.join(ROOT, 'inputs', 'alpha_isotype_alignment_TUBA1Bnumbering.csv')
TEMPLATE = os.path.join(ROOT, 'inputs', '9M1M.cif.gz')
GRAFT = os.path.join(ROOT, 'inputs', '4I4T.cif.gz')
AA3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS',
       'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL',
       'W': 'TRP', 'Y': 'TYR'}
LOOP = range(32, 46)


def forcefield():
    return app.ForceField(os.path.join(FF, 'amberff14SB.xml'), os.path.join(FF, 'gtp_meagher.xml'),
                          os.path.join(FF, 'water_ions_tip3p_iod.xml'))


def load_alignment():
    rows = list(csv.DictReader(open(ALN)))
    return {iso: {int(r['pos_TUBA1B']): r[iso] for r in rows} for iso in ('Tuba1a', 'Tuba1b', 'Tuba4a', 'Tuba8')}


def target_sequence(spec, aln):
    """Return {pos_TUBA1B: one-letter} for the construct, and a list of applied changes."""
    iso, _, mods = spec.partition(':')
    seq = {p: a for p, a in aln[iso].items() if a != '-'}
    changes = []
    for m in [x for x in mods.split(',') if x]:
        if m.startswith('loop='):
            donor = m.split('=')[1]
            for p in LOOP:
                if seq[p] != aln[donor][p]:
                    changes.append(f'{seq[p]}{p}{aln[donor][p]}')
                seq[p] = aln[donor][p]
        else:
            wt, pos, mt = m[0], int(m[1:-1]), m[-1]
            assert seq[pos] == wt, f'{spec}: position {pos} is {seq[pos]}, not {wt}'
            seq[pos] = mt
            changes.append(m)
    return iso, seq, changes


def octahedral_waters(mg, ligO, protein_O):
    """Positions (A) of water O atoms that complete an octahedron around mg given two bound ligand oxygens."""
    u1 = (ligO[0] - mg) / np.linalg.norm(ligO[0] - mg)
    u2 = ligO[1] - mg
    e2 = u2 - np.dot(u2, u1) * u1
    e2 /= np.linalg.norm(e2)
    e3 = np.cross(u1, e2)
    sites = [-u1, -e2, e3, -e3]
    out = []
    for s in sites:
        p = mg + 2.08 * s
        if len(protein_O) and np.min(np.linalg.norm(protein_O - p, axis=1)) < 2.4:
            continue  # a protein oxygen already sits there
        out.append(p)
    return out


def water_atoms(o, mg):
    """O, H1, H2 for a water whose dipole points away from Mg."""
    d = (o - mg) / np.linalg.norm(o - mg)
    a = np.cross(d, [1.0, 0, 0])
    if np.linalg.norm(a) < 1e-3:
        a = np.cross(d, [0, 1.0, 0])
    a /= np.linalg.norm(a)
    half = math.radians(104.52 / 2)
    h1 = o + 0.9572 * (math.cos(half) * d + math.sin(half) * a)
    h2 = o + 0.9572 * (math.cos(half) * d - math.sin(half) * a)
    return o, h1, h2


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1); v = b0 - b0 @ b1 * b1; w = b2 - b2 @ b1 * b1
    return float(np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w)))


def d_residues(top, X):
    """protein residues whose N-CA-C-CB dihedral is positive, i.e. D (about -124 deg for L)"""
    out = []
    for r in top.residues():
        A = {at.name: at.index for at in r.atoms()}
        if r.name != 'GLY' and all(k in A for k in ('N', 'CA', 'C', 'CB')) and _dihedral(X[A['N']], X[A['CA']], X[A['C']], X[A['CB']]) > 0:
            out.append(r)
    return out


def ideal_cb(n, ca, c):
    """ideal L CB from backbone N, CA, C (Angstrom)"""
    b = ca - n; cc = c - ca; aa = np.cross(b, cc)
    return -0.58273431 * aa + 0.56802827 * b - 0.54067466 * cc + ca


def enforce_L(fixer, rounds=3):
    """PDBFixer builds a side chain by superposing a template residue on the atoms that exist and then relaxing the
    new atoms; where only N, CA, C, O exist (a template Gly mutated to another residue, e.g. TUBA8 G43I / G45D on the
    TUBA1B backbone) CB can end up on the D side (2026-11-05: D-Ile43 / D-Asp45 in seven TUBA8-loop builds).  For each
    D residue: keep N, CA, C, O (and OXT), put CB at its ideal L position, and let PDBFixer rebuild the rest of the
    side chain on N/CA/C/O/CB.  Protein-only topologies (no hydrogens) as at this stage of a01 / a07.
    Returns (fixer, list of residues rebuilt)."""
    fixed = []
    for _ in range(rounds):
        X = np.array(fixer.positions.value_in_unit(unit.angstrom))
        bad = d_residues(fixer.topology, X)
        if not bad:
            return fixer, fixed
        badidx = {r.index for r in bad}
        fixed += [f'{r.chain.id}:{r.name}{r.id}' for r in bad]
        buf = io.StringIO()
        app.PDBFile.writeFile(fixer.topology, fixer.positions, buf, keepIds=True)
        atoms = list(fixer.topology.atoms()); k = 0; out = []
        for l in buf.getvalue().splitlines():
            if l.startswith(('ATOM', 'HETATM')):
                at = atoms[k]; k += 1; r = at.residue
                if r.index in badidx and at.name not in ('N', 'CA', 'C', 'O', 'OXT'):
                    continue
                out.append(l)
                if r.index in badidx and at.name == 'CA':
                    A = {x.name: x.index for x in r.atoms()}
                    cb = ideal_cb(X[A['N']], X[A['CA']], X[A['C']])
                    out.append(l[:12] + ' CB ' + l[16:30] + f'{cb[0]:8.3f}{cb[1]:8.3f}{cb[2]:8.3f}' + l[54:])
            else:
                out.append(l)
        assert k == len(atoms), (k, len(atoms))
        new = PDBFixer(pdbfile=io.StringIO('\n'.join(out) + '\n'))
        new.missingResidues = {}
        new.findMissingAtoms(); new.addMissingAtoms()
        fixer = new
    bad = d_residues(fixer.topology, np.array(fixer.positions.value_in_unit(unit.angstrom)))
    if bad:
        raise SystemExit(f'D residues remain after {rounds} rebuilds: {[f"{r.chain.id}:{r.name}{r.id}" for r in bad]}')
    return fixer, fixed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sid'); ap.add_argument('spec')
    ap.add_argument('--pad', type=float, default=1.0)
    ap.add_argument('--tail', choices=['none', 'full'], default='none',
                    help='none: stop at the last resolved residue (437), the disordered C-terminal tail is studied in the CAP-Gly-tail systems')
    ap.add_argument('--out', default=os.path.join(ROOT, 'build'))
    a = ap.parse_args()
    W = os.path.join(a.out, a.sid); os.makedirs(W, exist_ok=True)
    aln = load_alignment()
    iso, seq, changes = target_sequence(a.spec, aln)

    st = gemmi.read_structure(TEMPLATE); st.setup_entities()
    ch = st[0]['a']
    poly = ch.get_polymer()
    tmpl_seq = {r.seqid.num: gemmi.find_tabulated_residue(r.name).one_letter_code.upper() for r in poly}
    assert all(tmpl_seq[p] == aln['Tuba1b'][p] for p in tmpl_seq), 'template is not TUBA1B'
    last = max(tmpl_seq)
    # --- protein, threaded
    out = gemmi.Structure(); out.add_model(gemmi.Model('1')); oc = gemmi.Chain('A')
    nmut = 0
    for r in poly:
        p = r.seqid.num
        want = seq[p]
        nr = gemmi.Residue(); nr.seqid = gemmi.SeqId(p, ' '); nr.het_flag = 'A'
        keep = None
        if want != tmpl_seq[p]:
            nr.name = AA3[want]; nmut += 1
            keep = {'N', 'CA', 'C', 'O'} | ({'CB'} if 'G' not in (want, tmpl_seq[p]) else set())
        else:
            nr.name = r.name
        for at in r:
            if at.element.name == 'H' or (keep is not None and at.name not in keep):
                continue
            nr.add_atom(at)
        oc.add_residue(nr)
    out[0].add_chain(oc)
    out.setup_entities()
    pdb_txt = out.make_pdb_string()
    fixer = PDBFixer(pdbfile=io.StringIO(pdb_txt))
    tail = [AA3[seq[p]] for p in sorted(seq) if p > last] if a.tail == 'full' else []
    if tail:
        fixer.missingResidues = {(0, len(list(fixer.topology.residues()))): tail}
    else:
        fixer.missingResidues = {}
        seq = {p: x for p, x in seq.items() if p <= last}
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer, rebuilt_L = enforce_L(fixer)
    # --- GTP + Mg2+ + its four inner-shell waters, grafted from the 1.8 A crystal 4I4T chain A (canonical N-site:
    #     Mg bound to O1B, O1G and four waters).  Superposed on the CA atoms of the N-site residues.
    app.Topology.loadBondDefinitions(os.path.join(FF, 'gtp_bonds.xml'))
    ref = gemmi.read_structure(GRAFT); ref.setup_entities(); rch = ref[0]['A']
    rg = [r for r in rch if r.name == 'GTP'][0]; rmg = [r for r in rch if r.name == 'MG'][0]
    rmgpos = np.array(rmg[0].pos.tolist())
    rgp = np.array([at.pos.tolist() for at in rg if at.element.name != 'H'])
    site = sorted({r.seqid.num for r in rch.get_polymer() for at in r
                   if np.min(np.linalg.norm(rgp - np.array(at.pos.tolist()), axis=1)) < 6.0})
    tca = {r.seqid.num: np.array(r['CA'][0].pos.tolist()) for r in poly if r.find_atom('CA', '*')}
    rca = {r.seqid.num: np.array(r['CA'][0].pos.tolist()) for r in rch.get_polymer() if r.find_atom('CA', '*')}
    common = [p for p in site if p in tca and p in rca]
    X = np.array([rca[p] for p in common]); Y = np.array([tca[p] for p in common])
    xc, yc = X.mean(0), Y.mean(0)
    U, S_, Vt = np.linalg.svd((X - xc).T @ (Y - yc))
    dd = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, dd]) @ U.T
    move = lambda x: (R @ (np.asarray(x) - xc)) + yc
    site_rmsd = float(np.sqrt(np.mean(np.sum((np.array([move(x) for x in X]) - Y) ** 2, 1))))
    tg = [r for r in ch.get_ligands() if r.name == 'GTP'][0]
    tgp = {at.name: np.array(at.pos.tolist()) for at in tg}
    gtp_rmsd = float(np.sqrt(np.mean([np.sum((move(at.pos.tolist()) - tgp[at.name]) ** 2) for at in rg if at.name in tgp])))
    ns = gemmi.NeighborSearch(ref[0], ref.cell, 4).populate()
    shell = []
    for mk in ns.find_atoms(rmg[0].pos, '\0', radius=2.6):
        cra = mk.to_cra(ref[0])
        if cra.residue.name == 'HOH':
            shell.append(move(cra.atom.pos.tolist()))
    assert len(shell) == 4, f'expected 4 Mg-bound waters in 4I4T, got {len(shell)}'
    mgpos = move(rmgpos)
    lines = []
    k = 1
    for at in rg:
        if at.element.name == 'H':
            continue  # 4I4T carries riding hydrogens with CCD names; rebuilt below with the template's names
        x = move(at.pos.tolist())
        lines.append(f"HETATM{k:5d} {at.name:<4s} GTP B 501    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 20.00          {at.element.name:>2s}"); k += 1
    lines.append(f"HETATM{k:5d} MG    MG C 502    {mgpos[0]:8.3f}{mgpos[1]:8.3f}{mgpos[2]:8.3f}  1.00 20.00          MG"); k += 1
    for i, o in enumerate(shell):
        for nm, x in zip(('O', 'H1', 'H2'), water_atoms(o, mgpos)):
            lines.append(f"HETATM{k:5d} {nm:<4s} HOH D {601 + i:3d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00 20.00           {nm[0]}"); k += 1
    protH = np.array([np.array(at.pos.tolist()) for r in poly for at in r if at.element.name != 'H'])
    wsites = shell
    clash = float(min(np.min(np.linalg.norm(protH - o, axis=1)) for o in shell + [mgpos]))
    lig_pdb = app.PDBFile(io.StringIO('\n'.join(lines) + '\nEND\n'))
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.add(lig_pdb.topology, lig_pdb.positions)
    with open(os.path.join(W, 'start_noH.pdb'), 'w') as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    ff = forcefield()
    app.Modeller.loadHydrogenDefinitions(os.path.join(FF, 'gtp_hydrogens.xml'))
    modeller.addHydrogens(ff, pH=7.0)
    dres = d_residues(modeller.topology, np.array(modeller.positions.value_in_unit(unit.angstrom)))
    if dres:
        raise SystemExit(f'D residues after protonation: {[f"{r.chain.id}:{r.name}{r.id}" for r in dres]}')
    his = [r.name for r in modeller.topology.residues() if r.name in ('HIS', 'HID', 'HIE', 'HIP')]
    # --- solvate
    modeller.addSolvent(ff, model='tip3p', padding=a.pad * unit.nanometer, boxShape='dodecahedron',
                        ionicStrength=0.15 * unit.molar, positiveIon='Na+', negativeIon='Cl-', neutralize=True)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds)
    nbf = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
    qtot = sum(nbf.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))
    seqstr = ''.join(seq[p] for p in sorted(seq))
    with gzip.open(os.path.join(W, 'solvated.cif.gz'), 'wt') as fh:
        app.PDBxFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    resn = {}
    for r in modeller.topology.residues():
        resn[r.name] = resn.get(r.name, 0) + 1
    meta = dict(sid=a.sid, spec=a.spec, isoform=iso, changes=changes, n_rethreaded=nmut, length=len(seq),
                sequence=seqstr, tail_built=len(tail), n_atoms=system.getNumParticles(),
                net_charge=round(qtot, 6), box_nm=[[round(x, 4) for x in v.value_in_unit(unit.nanometer)]
                                                   for v in modeller.topology.getPeriodicBoxVectors()],
                mg_shell_waters_added=len(wsites), graft='4I4T chain A', graft_site_ca_rmsd_A=round(site_rmsd, 3), graft_vs_9M1M_gtp_rmsd_A=round(gtp_rmsd, 3), graft_min_dist_to_protein_heavy_A=round(clash, 2), tail_mode=a.tail, residue_counts={k: v for k, v in resn.items() if k in ('HOH', 'NA', 'CL', 'MG', 'GTP')},
                side_chains_rebuilt_L=rebuilt_L, histidines=his, forcefield=['amberff14SB.xml (ISOLDE)', 'gtp_meagher.xml', 'water_ions_tip3p_iod.xml'],
                template='9M1M chain a', openmm=mm.__version__)
    json.dump(meta, open(os.path.join(W, 'build.json'), 'w'), indent=1)
    print(json.dumps({k: meta[k] for k in ('sid', 'changes', 'length', 'n_atoms', 'net_charge', 'mg_shell_waters_added', 'graft_site_ca_rmsd_A', 'graft_vs_9M1M_gtp_rmsd_A', 'graft_min_dist_to_protein_heavy_A', 'residue_counts')}))


if __name__ == '__main__':
    main()
