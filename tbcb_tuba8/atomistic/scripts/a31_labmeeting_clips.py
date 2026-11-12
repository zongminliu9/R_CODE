"""Lab-meeting material from finished trajectories: thinned, superposed heavy-atom clips and per-frame series.

For each entry of the config (a JSON list) this writes, under --out:
  <label>.pdb / <label>.dcd   protein heavy atoms of every chain + GTP + Mg (no H, no water), every frame superposed on
                              the alpha core CA (alpha CA outside 25-55 and resSeq <= 437, the a26 core without the tail),
                              thinned to n_out frames between t0_ns and t1_ns (whole run by default)
  <label>_series.json         full-resolution series used for the figures:
     loop32_45  Q(lambda = 1.2) and backbone RMSD after core superposition, exactly as a26_loop_cvs.py
     res114     packing of residue 114: heavy atoms of other residues within 4.5 A of its side chain (CB and beyond);
                for SER also d(OG114 - O110), the hydrogen bond reported in the residue analysis
     groove     (kind 'tbcb') centroid distances from the tail aromatic ring (alpha PHE449 or TYR451) to the TBCB groove
                rings F187, W192 (six-membered ring), F223 and to K203 NZ; for TYR also d(OH451 - NZ K188)
     units      (kind 'md450', every replica r*) Q_u(t) for the folding units, the 32-45 loop and 110-120, with the
                a18_unfold_analysis.py definitions (contacts of the minimised start, |i-j| > 3, 4.5 A, beta 50/nm, 1.2)
A clip is written only for the replica named in the entry; mirror-image (D) alpha chains are refused, not reflected.
Usage: a31_labmeeting_clips.py <config.json> --out <dir>
Config entry: {"label": ..., "dir": <run dir with solute.pdb and seg_*.dcd>, "chain": "A", "kind": "md310|tbcb|md450",
               "n_out": 400, "t0_ns": null, "t1_ns": null, "clip": true}  (md450: "dir" is the system dir, "clip_rep": "r0";
               "clip": false writes the series only)
"""
import os, sys, json, glob, argparse
import numpy as np
import mdtraj as md

ap = argparse.ArgumentParser(); ap.add_argument('config'); ap.add_argument('--out', required=True)
ap.add_argument('--force', action='store_true', help='redo entries whose <label>_series.json exists')
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
FU = json.load(open(os.path.join(HERE, '..', 'inputs', 'folding_units.json')))
unit_of = {r: u for r, u in zip(FU['residues'], FU['unit']) if r <= 437}
UNITS = {f'unit{u}': sorted(r for r, x in unit_of.items() if x == u) for u in sorted(set(unit_of.values()))}
UNITS['loop32_45'] = list(range(32, 46)); UNITS['seg110_120'] = list(range(110, 121))


def chains_of(pdb):
    return np.array([l[21] for l in open(pdb) if l.startswith(('ATOM', 'HETATM'))])


def load(rdir):
    pdb = os.path.join(rdir, 'solute.pdb')
    segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd')))
    assert segs, f'no trajectory in {rdir}'
    t = md.join([md.load_dcd(s, top=pdb) for s in segs], check_topology=False)
    ref = md.load(pdb)
    ns = None
    for f in ('progress.json', 'status.json'):
        p = os.path.join(rdir, f)
        if os.path.exists(p):
            ns = json.load(open(p)).get('ns_done'); break
    assert ns, f'no ns_done in {rdir}'
    return t, ref, pdb, float(ns) / t.n_frames, segs


def alpha_atoms(top, ch, chain):
    return [x for x in top.atoms if x.residue.is_protein and ch[x.index] == chain]


def handedness(ref, alpha):
    by = {(x.residue.resSeq, x.name): x.index for x in alpha}
    res = sorted({x.residue.resSeq for x in alpha if x.residue.name != 'GLY'})
    quad = np.array([[by[(r, n)] for n in ('N', 'CA', 'C', 'CB')] for r in res if all((r, n) in by for n in ('N', 'CA', 'C', 'CB'))])
    return float(np.mean(md.compute_dihedrals(ref, quad)[0] < 0))


def loop_series(t, ref, alpha):
    rid = {x.index: x.residue.resSeq for x in alpha}
    loop = [x.index for x in alpha if 32 <= x.residue.resSeq <= 45 and x.element.symbol != 'H']
    rest = [x.index for x in alpha if not (32 <= x.residue.resSeq <= 45) and x.element.symbol != 'H']
    P0 = ref.xyz[0]
    d = np.linalg.norm(P0[loop][:, None] - P0[rest][None], axis=2)
    ii, jj = np.where(d < 0.45)
    keep = [k for k in range(len(ii)) if abs(rid[loop[ii[k]]] - rid[rest[jj[k]]]) > 3]
    pairs = np.array([[loop[ii[k]], rest[jj[k]]] for k in keep]); d0 = d[ii[keep], jj[keep]]
    D = md.compute_distances(t, pairs, periodic=False)
    Q = (1.0 / (1.0 + np.exp(50.0 * (D - 1.2 * d0)))).mean(axis=1)
    core = [x.index for x in alpha if x.name == 'CA' and not (25 <= x.residue.resSeq <= 55) and x.residue.resSeq <= 437]
    lbb = [x.index for x in alpha if x.name in ('N', 'CA', 'C', 'O') and 32 <= x.residue.resSeq <= 45]
    ts = md.Trajectory(t.xyz.copy(), t.topology); ts.superpose(ref, 0, atom_indices=core)
    rmsd = np.sqrt(np.mean(np.sum((ts.xyz[:, lbb] - P0[lbb]) ** 2, axis=2), axis=1)) * 10
    return {'Q_lam1.2': np.round(Q, 3).tolist(), 'rmsd_loop_A': np.round(rmsd, 2).tolist(), 'n_native_pairs': int(len(pairs))}, ts, core


def res114_series(t, ref, alpha):
    by = {(x.residue.resSeq, x.name): x.index for x in alpha}
    r114 = [x for x in alpha if x.residue.resSeq == 114][0].residue
    sc = [x.index for x in alpha if x.residue.resSeq == 114 and x.element.symbol != 'H' and x.name not in ('N', 'CA', 'C', 'O')]
    own = set(x.index for x in alpha if x.residue.resSeq == 114)
    heavy = np.array([x.index for x in t.topology.atoms if x.element is not None and x.element.symbol != 'H'
                      and x.residue.name not in ('HOH', 'WAT') and x.index not in own])
    P0 = ref.xyz[0]
    near = heavy[np.min(np.linalg.norm(P0[heavy][:, None] - P0[sc][None], axis=2), axis=1) < 1.2]
    pairs = np.array([[i, j] for i in sc for j in near])
    D = md.compute_distances(t, pairs, periodic=False).reshape(t.n_frames, len(sc), len(near))
    ncont = (D.min(axis=1) < 0.45).sum(axis=1)
    out = {'resname': r114.name, 'n_contacts_sc_4.5A': ncont.astype(int).tolist()}
    if r114.name == 'SER':
        out['d_OG114_O110_A'] = np.round(md.compute_distances(t, [[by[(114, 'OG')], by[(110, 'O')]]], periodic=False)[:, 0] * 10, 2).tolist()
    return out


RING6 = {'PHE': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'), 'TYR': ('CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'),
         'TRP': ('CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2')}


def groove_series(t, ch, chain_a, chain_b):
    top = t.topology
    resmap = {}
    for x in top.atoms:
        if x.residue.is_protein:
            resmap.setdefault((ch[x.index], x.residue.resSeq), (x.residue.name, {}))[1][x.name] = x.index
    atoms = lambda chain, resseq: resmap[(chain, resseq)]
    tail = next(((rs,) + resmap[(chain_a, rs)] for rs in (451, 449)
                 if (chain_a, rs) in resmap and resmap[(chain_a, rs)][0] in ('PHE', 'TYR')), None)
    assert tail, 'no aromatic tail residue at 449/451'
    cen = lambda ix, names: t.xyz[:, [ix[n] for n in names]].mean(axis=1)
    c_tail = cen(tail[2], RING6[tail[1]])
    out = {'tail_residue': f'{tail[1]}{tail[0]}'}
    for rs in (187, 192, 223):
        nm, ix = atoms(chain_b, rs)
        out[f'd_ring_{nm}{rs}_A'] = np.round(np.linalg.norm(c_tail - cen(ix, RING6[nm]), axis=1) * 10, 2).tolist()
    nm, ix = atoms(chain_b, 203)
    out[f'd_ring_NZ_{nm}203_A'] = np.round(np.linalg.norm(c_tail - t.xyz[:, ix['NZ']], axis=1) * 10, 2).tolist()
    if tail[1] == 'TYR':
        nm, ix = atoms(chain_b, 188)
        if 'NZ' in ix:
            out['d_OH451_NZ188_A'] = np.round(np.linalg.norm(t.xyz[:, tail[2]['OH']] - t.xyz[:, ix['NZ']], axis=1) * 10, 2).tolist()
    return out


def unit_series(t, ref, top):
    heavy = top.select('protein and not element H')
    rs = np.array([top.atom(i).residue.resSeq for i in heavy])
    X0 = ref.xyz[0][heavy]
    out = {}
    for g, res in UNITS.items():
        gi = np.where(np.isin(rs, res))[0]
        dd = np.linalg.norm(X0[gi][:, None] - X0[None], axis=2)
        ii, jj = np.where((dd < 0.45) & (np.abs(rs[gi][:, None] - rs[None]) > 3))
        pairs = np.stack([heavy[gi][ii], heavy[jj]], axis=1); r0 = dd[ii, jj]
        D = md.compute_distances(t, pairs, periodic=False)
        out[g] = np.round(np.mean(1.0 / (1.0 + np.exp(50.0 * (D - 1.2 * r0[None]))), axis=1), 3).tolist()
    return out


def write_clip(ts, label, e, ns_per_frame):
    top = ts.topology
    keep = [x.index for x in top.atoms if x.element is not None and x.element.symbol != 'H' and x.residue.name not in ('HOH', 'WAT', 'NA', 'CL', 'K')]
    n = ts.n_frames
    i0 = 0 if e.get('t0_ns') is None else int(round(e['t0_ns'] / ns_per_frame))
    i1 = n - 1 if e.get('t1_ns') is None else min(n - 1, int(round(e['t1_ns'] / ns_per_frame)))
    idx = np.unique(np.linspace(i0, i1, e.get('n_out', 400)).round().astype(int))
    clip = ts[idx].atom_slice(keep)
    clip[0].save_pdb(os.path.join(a.out, f'{label}.pdb'))
    clip.save_dcd(os.path.join(a.out, f'{label}.dcd'))
    return {'frames_full': int(n), 'clip_frame_index': idx.tolist(), 'clip_time_ns': np.round(idx * ns_per_frame, 3).tolist(),
            'clip_atoms': int(clip.n_atoms)}


summary = {}
def run_entry(e):
    label = e['label']
    rdirs = sorted(glob.glob(os.path.join(e['dir'], 'r*'))) if e['kind'] == 'md450' else [e['dir']]
    rec = {'entry': e, 'replicas': {}}
    for rdir in rdirs:
        rep = os.path.basename(rdir)
        t, ref, pdb, nspf, segs = load(rdir)
        ch = chains_of(pdb); assert len(ch) == t.n_atoms == ref.n_atoms
        alpha = alpha_atoms(t.topology, ch, e['chain'])
        fl = handedness(ref, alpha)
        r = {'ns_per_frame': nspf, 'n_frames': int(t.n_frames), 'segments': [os.path.basename(s) for s in segs], 'alpha_fraction_L': fl}
        if fl < 0.95:
            print('  refused: alpha chain not L', rep, fl, flush=True)
            r['refused'] = 'alpha chain not L'; rec['replicas'][rep] = r; continue
        r['loop32_45'], ts, core = loop_series(t, ref, alpha)
        r['res114'] = res114_series(t, ref, alpha)
        if e['kind'] == 'tbcb':
            r['groove'] = groove_series(t, ch, e['chain'], e['partner'])
        if e['kind'] == 'md450':
            r['units'] = unit_series(t, ref, t.topology)
        if e.get('clip', True) and (e['kind'] != 'md450' or rep == e.get('clip_rep', 'r0')):
            r['clip'] = write_clip(ts, label if e['kind'] != 'md450' else f'{label}_{rep}', e, nspf)
        rec['replicas'][rep] = r
        print('  ', rep, 'frames', t.n_frames, 'ns/frame', round(nspf, 4), 'L', round(fl, 3),
              'Q mean', round(float(np.mean(r['loop32_45']['Q_lam1.2'])), 3), flush=True)
        del t, ts
    json.dump(rec, open(os.path.join(a.out, f'{label}_series.json'), 'w'))
    return {k: {kk: vv for kk, vv in v.items() if kk in ('n_frames', 'ns_per_frame', 'alpha_fraction_L', 'refused')} for k, v in rec['replicas'].items()}


# each entry is written as soon as it is done; a failing entry is recorded and the next one runs
for e in json.load(open(a.config)):
    label = e['label']
    if not a.force and os.path.exists(os.path.join(a.out, f'{label}_series.json')):
        summary[label] = 'exists, skipped'; continue
    print('==', label, e['dir'], flush=True)
    try:
        summary[label] = run_entry(e)
    except Exception as ex:
        import traceback; traceback.print_exc()
        summary[label] = {'failed': repr(ex)}
json.dump(summary, open(os.path.join(a.out, 'clips_summary.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
