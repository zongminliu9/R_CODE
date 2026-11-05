"""Interface stability and inter-chain contacts along the unbiased 310 K MD of a complex.

Per replica, frames after --skip_ns (default 0, so drift from the start is visible), every --stride-th frame:
  * heavy-atom contacts (< 4.5 A, the definition used in a02_minimize.py) between partner A and partner B;
  * Q_if: fraction of the reference residue-residue interface contacts still present, the reference being the
    residue pairs in contact (any heavy-atom pair < 4.5 A) in >= 50 % of the frames of the first 1 ns;
  * backbone RMSD of each partner after superposing on that partner, and the CA RMSD of the reference interface
    residues after superposing on partner A (how far B moves on A);
  * inter-chain H-bonds (Baker-Hubbard) and salt bridges (carboxylate O to Lys NZ / Arg NE,NH1,NH2 < 4 A) with
    occupancy >= 0.30 over the second half of the window;
  * for each --focus residue (chain:resSeq), every residue of the other partner it touches (heavy atoms < 4.5 A)
    with occupancy >= 0.20, in each half of the window.
The saved coordinates are not wrapped (a03_md.py writes enforcePeriodicBox=False), so direct distances are the
physical ones; the minimum-image contact count is also computed and its largest difference is reported as a check.
Values are given for 4 contiguous time blocks so drift is visible.  These are structural observables, not free
energies.
Usage: a25_interface_md.py <md_root>/<sid> <chainsA> <chainsB> <out.json> [--focus A:188,A:203] [--stride 5]
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md

ap = argparse.ArgumentParser()
ap.add_argument('mdsys'); ap.add_argument('ca'); ap.add_argument('cb'); ap.add_argument('out')
ap.add_argument('--focus', default=''); ap.add_argument('--stride', type=int, default=5)
ap.add_argument('--skip_ns', type=float, default=0.0); ap.add_argument('--dt_ps', type=float, default=20.0)
a = ap.parse_args()
CHA, CHB = set(a.ca.split(',')), set(a.cb.split(','))
ACID = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}
BASE = {'LYS': ('NZ',), 'ARG': ('NE', 'NH1', 'NH2')}


def pdb_chains(fn):
    """chain letter of every atom, in file order (mdtraj keeps the PDB atom order)"""
    out = []
    for l in open(fn):
        if l.startswith(('ATOM', 'HETATM')):
            out.append(l[21])
    return np.array(out)


def blocks(x, nb=4):
    return [float(np.mean(b)) for b in np.array_split(np.asarray(x, float), nb)]


def rlabel(r, ch):
    return f'{ch}:{r.name}{r.resSeq}'


reps_out = {}
for rdir in sorted(glob.glob(os.path.join(a.mdsys, 'r*'))):
    segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd')))
    pdb = os.path.join(rdir, 'solute.pdb')
    if not segs or not os.path.exists(pdb):
        continue
    t = md.join([md.load_dcd(s, top=pdb) for s in segs], check_topology=False)
    t_end_ns = t.n_frames * a.dt_ps / 1000.0
    t = t[int(a.skip_ns * 1000 / a.dt_ps)::a.stride]
    if t.n_frames < 20:
        continue
    top = t.topology
    ch = pdb_chains(pdb)
    assert len(ch) == top.n_atoms, (len(ch), top.n_atoms)
    res_ch = {r.index: ch[next(iter(r.atoms)).index] for r in top.residues}
    prot = lambda chs: [x.index for x in top.atoms if x.residue.is_protein and ch[x.index] in chs]
    heavy = lambda idx: [i for i in idx if top.atom(i).element.symbol != 'H']
    A, B = prot(CHA), prot(CHB)
    HA, HB = np.array(heavy(A)), np.array(heavy(B))
    assert len(HA) and len(HB), 'empty partner selection'
    time_ns = (a.skip_ns + np.arange(t.n_frames) * a.stride * a.dt_ps / 1000.0)
    # candidate heavy-atom pairs: within 1.2 nm in the first, middle or last frame
    cand = set()
    for f in (0, t.n_frames // 2, t.n_frames - 1):
        XA, XB = t.xyz[f, HA], t.xyz[f, HB]
        for k in range(0, len(HA), 1500):
            d = np.linalg.norm(XA[k:k + 1500, None] - XB[None], axis=2)
            ii, jj = np.where(d < 1.2)
            cand |= set(zip(HA[k + ii].tolist(), HB[jj].tolist()))
    pairs = np.array(sorted(cand))
    D = md.compute_distances(t, pairs, periodic=False)
    Dp = md.compute_distances(t, pairs, periodic=True)
    n_contacts = (D < 0.45).sum(axis=1)
    n_contacts_mi = (Dp < 0.45).sum(axis=1)
    # residue-pair minimum distances
    ra = np.array([top.atom(i).residue.index for i in pairs[:, 0]]); rb = np.array([top.atom(j).residue.index for j in pairs[:, 1]])
    key = ra * top.n_residues + rb
    order = np.argsort(key, kind='stable'); key_s = key[order]
    starts = np.r_[0, np.where(np.diff(key_s))[0] + 1]
    ukeys = key_s[starts]
    Dres = np.minimum.reduceat(D[:, order], starts, axis=1)            # (frames, residue pairs)
    C = Dres < 0.45
    n0 = max(1, int(round(1.0 / (a.stride * a.dt_ps / 1000.0))))     # frames in the first 1 ns
    ref = C[:n0].mean(axis=0) >= 0.5
    Q = C[:, ref].mean(axis=1) if ref.any() else np.zeros(t.n_frames)
    # RMSDs
    bbA = [i for i in A if top.atom(i).name in ('N', 'CA', 'C')]; bbB = [i for i in B if top.atom(i).name in ('N', 'CA', 'C')]
    rmsdA = md.rmsd(t, t, 0, atom_indices=bbA) * 10; rmsdB = md.rmsd(t, t, 0, atom_indices=bbB) * 10
    ref_res = set((ukeys[ref] // top.n_residues).tolist()) | set((ukeys[ref] % top.n_residues).tolist())
    iface_ca = [x.index for x in top.atoms if x.name == 'CA' and x.residue.index in ref_res]
    caA = [i for i in A if top.atom(i).name == 'CA']
    ts = md.Trajectory(t.xyz.copy(), top); ts.superpose(ts, 0, atom_indices=caA)
    d_if = np.sqrt(np.mean(np.sum((ts.xyz[:, iface_ca] - ts.xyz[0, iface_ca]) ** 2, axis=2), axis=1)) * 10 if iface_ca else np.zeros(t.n_frames)
    half = t.n_frames // 2
    # inter-chain H-bonds, second half
    th = t[half:]
    hb = md.baker_hubbard(th, freq=0.30, periodic=False)
    hbl = []
    if len(hb):
        dist = md.compute_distances(th, hb[:, [1, 2]], periodic=False)
        ang = md.compute_angles(th, hb, periodic=False)
        occ = ((dist < 0.25) & (ang > 2.0 * np.pi / 3)).mean(axis=0)
        for (dd, hh, aa), o in zip(hb, occ):
            rd, ra_ = top.atom(dd).residue, top.atom(aa).residue
            if (res_ch[rd.index] in CHA) != (res_ch[ra_.index] in CHA) and rd.is_protein and ra_.is_protein:
                hbl.append([f'{rlabel(rd, res_ch[rd.index])}:{top.atom(dd).name}', f'{rlabel(ra_, res_ch[ra_.index])}:{top.atom(aa).name}', round(float(o), 3)])
    # inter-chain salt bridges, second half
    sbl = []
    resA = [r for r in top.residues if r.is_protein and res_ch[r.index] in CHA]
    resB = [r for r in top.residues if r.is_protein and res_ch[r.index] in CHB]
    for r1 in resA:
        for r2 in resB:
            if not ((r1.name in ACID and r2.name in BASE) or (r1.name in BASE and r2.name in ACID)):
                continue
            g1 = [x.index for x in r1.atoms if x.name in ACID.get(r1.name, ()) + BASE.get(r1.name, ())]
            g2 = [x.index for x in r2.atoms if x.name in ACID.get(r2.name, ()) + BASE.get(r2.name, ())]
            # quick screen on the first frame of the half
            if np.min(np.linalg.norm(th.xyz[0, g1][:, None] - th.xyz[0, g2][None], axis=2)) > 1.5:
                continue
            dmin = md.compute_distances(th, np.array([[i, j] for i in g1 for j in g2]), periodic=False).min(axis=1)
            o = float(np.mean(dmin < 0.40))
            if o >= 0.30:
                sbl.append([rlabel(r1, res_ch[r1.index]), rlabel(r2, res_ch[r2.index]), round(o, 3)])
    # focus residues
    foc = {}
    for item in [x for x in a.focus.split(',') if x]:
        c0, n = item.split(':'); n = int(n)
        rr = [r for r in top.residues if r.is_protein and res_ch[r.index] == c0 and r.resSeq == n]
        if not rr:
            foc[item] = 'not found'; continue
        r = rr[0]
        other = CHB if c0 in CHA else CHA
        mine = [x.index for x in r.atoms if x.element.symbol != 'H']
        partners = [q for q in top.residues if q.is_protein and res_ch[q.index] in other]
        rec = {'resname': r.name, 'first_half': [], 'second_half': []}
        for q in partners:
            th_ = [x.index for x in q.atoms if x.element.symbol != 'H']
            if np.min(np.linalg.norm(t.xyz[::max(1, t.n_frames // 20)][:, mine][:, :, None] - t.xyz[::max(1, t.n_frames // 20)][:, th_][:, None], axis=3)) > 1.0:
                continue
            dm = md.compute_distances(t, np.array([[i, j] for i in mine for j in th_]), periodic=False).min(axis=1)
            o1, o2 = float(np.mean(dm[:half] < 0.45)), float(np.mean(dm[half:] < 0.45))
            if o1 >= 0.20:
                rec['first_half'].append([rlabel(q, res_ch[q.index]), round(o1, 3)])
            if o2 >= 0.20:
                rec['second_half'].append([rlabel(q, res_ch[q.index]), round(o2, 3)])
        foc[item] = rec
    reps_out[os.path.basename(rdir)] = dict(
        n_frames=int(t.n_frames), trajectory_ns=float(t_end_ns), window_ns=[float(time_ns[0]), float(time_ns[-1])],
        block_edges_ns=[float(b[0]) for b in np.array_split(time_ns, 4)] + [float(time_ns[-1])],
        heavy_contacts_blocks=blocks(n_contacts), heavy_contacts_first_frame=int(n_contacts[0]),
        minimum_image_max_difference=int(np.max(np.abs(n_contacts_mi - n_contacts))),
        n_reference_residue_contacts=int(ref.sum()), Q_if_blocks=blocks(Q), Q_if_min=float(Q.min()),
        rmsd_bb_A_blocks_A=blocks(rmsdA), rmsd_bb_B_blocks_A=blocks(rmsdB), iface_ca_rmsd_on_A_blocks_A=blocks(d_if),
        interchain_hbonds_second_half=sorted(hbl, key=lambda x: -x[2]), interchain_salt_bridges_second_half=sorted(sbl, key=lambda x: -x[2]),
        focus=foc)
    print(os.path.basename(rdir), t.n_frames, 'frames', 'Q_if blocks', [round(x, 3) for x in reps_out[os.path.basename(rdir)]['Q_if_blocks']], flush=True)
json.dump({'system': os.path.basename(a.mdsys.rstrip('/')), 'chains_A': sorted(CHA), 'chains_B': sorted(CHB), 'replicas': reps_out},
          open(a.out, 'w'), indent=1)
print('wrote', a.out)
