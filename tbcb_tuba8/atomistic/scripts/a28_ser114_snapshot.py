"""Hydration and H-bond partners of the Ser114 hydroxyl in the last HREX configuration of every folded replicate.

I114S (TUBA1A): Ser is the appended ("new") side chain; its OG is identified among the appended atoms as the atom
  ~1.43 A from the appended CB with a partner ~0.96 A away (HG); states with Ser fully interacting: 19, 20, 21.
S114I (TUBA8): Ser is the WT residue 114 (atoms in wt_system.cif.gz); Ser fully interacting in states 0, 1, 2.
Counts: water O within 3.5 A of OG; protein N/O (other than Ser114 itself) within 3.5 A of OG; water O within 6 A
of the side-chain centroid (pocket water).  Minimum image with the equilibrated box (HREX runs at constant volume).
One configuration per replica: a snapshot, not an ensemble average.
Usage: a28_ser114_snapshot.py <dir holding a1a_114SER_folded/r*/ and a8_114ILE_folded/r*/ with fep_ck.npz,
       equil_state.npz and r0/wt_system.cif.gz copied from runs/fep>
"""
import gzip, json, os, sys
import numpy as np

D = sys.argv[1]


def read_cif(path):
    cols, rows, inloop = [], [], False
    for l in gzip.open(path, 'rt'):
        if l.startswith('_atom_site.'):
            cols.append(l.strip().split('.', 1)[1]); inloop = True; continue
        if inloop:
            if l.startswith(('ATOM', 'HETATM')):
                rows.append(l.split())
            elif cols and rows and (l.startswith('#') or l.startswith('loop_') or l.startswith('_')):
                break
    ix = {c: i for i, c in enumerate(cols)}
    return [dict(el=r[ix['type_symbol']], name=r[ix['label_atom_id']], res=r[ix['label_comp_id']],
                 seq=r[ix.get('auth_seq_id', ix['label_seq_id'])], chain=r[ix.get('auth_asym_id', ix['label_asym_id'])])
            for r in rows]


SHIFTS = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)


def mic(d, box):
    """Minimum image for OpenMM's reduced triclinic box (rows a, b, c): reduce along c, b, a, then take the shortest
    of the 27 neighbouring images (exact for any distance shorter than half the smallest box width)."""
    box = np.asarray(box, float)
    d = np.atleast_2d(d).copy()
    d -= np.outer(np.round(d[:, 2] / box[2, 2]), box[2])
    d -= np.outer(np.round(d[:, 1] / box[1, 1]), box[1])
    d -= np.outer(np.round(d[:, 0] / box[0, 0]), box[0])
    cand = d[:, None, :] + (SHIFTS @ box)[None]
    k = np.argmin((cand ** 2).sum(-1), axis=1)
    return cand[np.arange(len(d)), k]


for leg, ser_states, ser_is_new in (('a1a_114SER_folded', (19, 20, 21), True), ('a8_114ILE_folded', (0, 1, 2), False)):
    atoms = read_cif(os.path.join(D, leg, 'r0', 'wt_system.cif.gz'))
    nwt = len(atoms)
    wat_O = np.array([i for i, a in enumerate(atoms) if a['res'] in ('HOH', 'WAT') and a['el'] == 'O'])
    prot_NO = np.array([i for i, a in enumerate(atoms) if a['res'] not in ('HOH', 'WAT', 'NA', 'CL', 'K', 'MG', 'Na+', 'Cl-')
                        and a['el'] in ('N', 'O') and not (a['seq'] == '114' and a['res'] in ('SER', 'ILE'))])
    res114 = [i for i, a in enumerate(atoms) if a['seq'] == '114' and a['res'] in ('SER', 'ILE')]
    print(f'== {leg}: {nwt} WT atoms, {len(wat_O)} waters; residue 114 in WT topology:',
          {atoms[i]['name']: i for i in res114}.keys().__len__(), atoms[res114[0]]['res'])
    for rep in ('r0', 'r1', 'r2'):
        z = np.load(os.path.join(D, leg, rep, 'fep_ck.npz'), allow_pickle=True)
        pos, s, it = z['pos'], z['s'], int(z['iteration'])
        box = np.load(os.path.join(D, leg, rep, 'equil_state.npz'))['box']
        nnew = pos.shape[1] - nwt
        if ser_is_new:
            app = np.arange(nwt, nwt + nnew)
            X0 = pos[0]
            # identify OG among appended atoms: O bonded to CB (~0.143 nm) and HG (~0.096 nm)
            dd = np.linalg.norm(X0[app][:, None] - X0[app][None], axis=2)
            og = None
            for k in range(nnew):
                near = sorted(dd[k][dd[k] > 0])
                if len(near) >= 2 and abs(near[0] - 0.096) < 0.01 and abs(near[1] - 0.143) < 0.01:
                    og = app[k]
            sc = app
        else:
            og = [i for i in res114 if atoms[i]['name'] == 'OG'][0]
            sc = np.array([i for i in res114 if atoms[i]['name'] in ('CB', 'OG')])
        out = []
        for st in ser_states:
            r = int(np.where(s == st)[0][0]); X = pos[r]
            dW = np.linalg.norm(mic(X[wat_O] - X[og], box), axis=1)
            dP = np.linalg.norm(mic(X[prot_NO] - X[og], box), axis=1)
            cen = X[sc].mean(0)
            dWc = np.linalg.norm(mic(X[wat_O] - cen, box), axis=1)
            partners = [f"{atoms[prot_NO[k]]['res']}{atoms[prot_NO[k]]['seq']}:{atoms[prot_NO[k]]['name']} {dP[k]*10:.2f}"
                        for k in np.argsort(dP)[:3] if dP[k] < 0.35]
            out.append(f"state {st}: water O <3.5 A of OG {int((dW < 0.35).sum())}, waters <6 A of side chain {int((dWc < 0.6).sum())}, "
                       f"nearest water {dW.min()*10:.2f} A; protein N/O <3.5 A: {partners}")
        print(f'  {rep} (iteration {it}, {nnew} appended atoms, OG index {og})')
        for o in out:
            print('     ' + o)
