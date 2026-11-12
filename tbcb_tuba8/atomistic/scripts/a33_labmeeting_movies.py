"""Lab-meeting trajectories for PyMOL, MP4 movies and stills, from the a31 clips.

1. every clip is superposed on one common reference (alpha core CA, resSeq outside 25-55 and <= 437, matched by
   residue number; reference = frame 0 of free TUBA1A r1) so that side-by-side panels share one camera, and written
   as <out>/<section>/pymol/<name>.pdb + .dcd; a styled PyMOL session <name>.pse holds the same trajectory.
2. movies: PyMOL ray-traced frames (light 5-frame smoothing for display only; the .dcd/.pse keep the raw frames),
   composed with a title, a clock and a running read-out of the quantity the movie is about, encoded as H.264 MP4.
3. stills: representative frames picked from the a31 series (rules printed with each still).
Usage: a33_labmeeting_movies.py --out <labmeeting dir> --work <labmeeting_work dir> [--only movie1,movie2] [--no_movies]
"""
import os, sys, json, glob, argparse, subprocess, shutil
import numpy as np
import mdtraj as md
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--work', required=True)
ap.add_argument('--only', default=''); ap.add_argument('--no_movies', action='store_true'); ap.add_argument('--stills_only', action='store_true'); ap.add_argument('--with_stills', action='store_true'); ap.add_argument('--skip_prep', action='store_true'); ap.add_argument('--prep_names', default='')
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__)); SCENE = os.environ.get('A33_SCENE', os.path.join(HERE, 'a33_pymol_scene.py'))
PYMOL = shutil.which('pymol'); FFMPEG = shutil.which('ffmpeg')
FONT = '/System/Library/Fonts/Supplemental/Arial.ttf'; FONTB = '/System/Library/Fonts/Supplemental/Arial Bold.ttf'
BLUE, ORANGE, GRAY, INK, INK2, GRID = '#2a78d6', '#eb6834', '#898781', '#0b0b0b', '#52514e', '#e1e0d9'
ISOC = {'tuba1a': BLUE, 'tuba8': ORANGE}
TH = json.load(open(os.path.join(HERE, '..', 'interim', '20261104_2110', 'rest2_state_thresholds.json')))
QN, RN, QU = TH['Q12_p05'], TH['rmsd_p95_A'], TH['U_threshold_Q12']
FU = json.load(open(os.path.join(HERE, '..', 'inputs', 'folding_units.json')))
unit_of = {r: u for r, u in zip(FU['residues'], FU['unit']) if r <= 437}
UNITS = [sorted(r for r, x in unit_of.items() if x == u) for u in sorted(set(unit_of.values()))]

# name: (section, source dir, alpha chain, partner chain, isoform colour, style)
CL = {'free_TUBA1A_r1': ('04_loop32-45_states', 'out_marlowe', 'A', None, 'tuba1a', 'loop'),
      'free_TUBA1A_r0': ('04_loop32-45_states', 'out_sherlock', 'A', None, 'tuba1a', 'loop'),
      'free_TUBA8_r0': ('04_loop32-45_states', 'out_sherlock', 'A', None, 'tuba8', 'loop'),
      'heterodimer_TUBA1A_r1': ('04_loop32-45_states', 'out_sherlock', 'A', 'B', 'tuba1a', 'loop'),
      'heterodimer_TUBA8_r0': ('04_loop32-45_states', 'out_sherlock', 'A', 'B', 'tuba8', 'loop'),
      'TBCE9M1L_TUBA1A_r0': ('04_loop32-45_states', 'out_sherlock', 'B', 'A', 'tuba1a', 'loop'),
      'TBCE9M1L_TUBA8_r1': ('04_loop32-45_states', 'out_sherlock', 'B', 'A', 'tuba8', 'loop'),
      'TBCB_TUBA8_r1': ('02_TBCB_tail_clamp', 'out_sherlock', 'A', 'B', 'tuba8', 'groove'),
      'TBCB_TUBA1A_r2': ('02_TBCB_tail_clamp', 'out_sherlock', 'A', 'B', 'tuba1a', 'groove'),
      'md450_TUBA8_r0': ('06_450K_unfolding_order', 'out_sherlock', 'A', None, 'tuba8', 'unfold'),
      'md450_TUBA1A_r0': ('06_450K_unfolding_order', 'out_sherlock', 'A', None, 'tuba1a', 'unfold')}
W_ = os.path.join(a.work, 'render'); os.makedirs(W_, exist_ok=True)


def src(name, ext):
    return os.path.join(a.work, CL[name][1], f'{name}.{ext}')


def series(name):
    lab = name[:-3] if name.startswith('md450_') else name
    d = json.load(open(os.path.join(a.work, CL[name][1], f'{lab}_series.json')))
    rep = name[-2:] if name.startswith('md450_') else list(d['replicas'])[0]
    return d['replicas'][rep]


def core_idx(top, chain_letter_index):
    ch = top.chain(chain_letter_index)
    return {a_.residue.resSeq: a_.index for a_ in ch.atoms if a_.name == 'CA' and not (25 <= a_.residue.resSeq <= 55) and a_.residue.resSeq <= 437}


def chain_index(pdb, letter):
    seen = []
    for l in open(pdb):
        if l.startswith(('ATOM', 'HETATM')) and l[21] not in seen:
            seen.append(l[21])
    return seen.index(letter)


def superpose_all():
    ref_name = 'free_TUBA1A_r1'; rp = src(ref_name, 'pdb'); ref = md.load(rp)
    rc = core_idx(ref.topology, chain_index(rp, 'A'))
    done = []
    for name, (sec, _, A, B, iso, sty) in CL.items():
        if not os.path.exists(src(name, 'dcd')):
            continue
        d = os.path.join(a.out, sec, 'pymol'); os.makedirs(d, exist_ok=True)
        if a.prep_names and name not in a.prep_names.split(','):
            continue
        t = md.load_dcd(src(name, 'dcd'), top=src(name, 'pdb'))
        if sty == 'groove':
            # TBCB slides on alpha (interface CA RMSD up to ~14 A), so a groove clip is superposed on the TBCB CAP-Gly
            # domain (chain B residues 160-244, CA) of its own first frame: the groove stays put, alpha moves around it
            capg = [a_.index for a_ in t.topology.chain(chain_index(src(name, 'pdb'), B)).atoms if a_.name == 'CA' and 160 <= a_.residue.resSeq <= 244]
            t.superpose(t, 0, atom_indices=capg)
            rms = md.rmsd(t, t, 0, atom_indices=capg) * 10
            t[0].save_pdb(os.path.join(d, f'{name}.pdb')); t.save_dcd(os.path.join(d, f'{name}.dcd'))
            done.append(name); print(f'  superposed {name} on the TBCB CAP-Gly domain: {t.n_frames} frames, CAP-Gly CA RMSD to frame 0 median {np.median(rms):.2f} A', flush=True)
            continue
        mc = core_idx(t.topology, chain_index(src(name, 'pdb'), A))
        common = sorted(set(rc) & set(mc)); assert len(common) > 300, (name, len(common))
        t.superpose(ref, 0, atom_indices=[mc[r] for r in common], ref_atom_indices=[rc[r] for r in common])
        rms = md.rmsd(t[0], ref, 0, atom_indices=[mc[r] for r in common], ref_atom_indices=[rc[r] for r in common])[0] * 10
        t[0].save_pdb(os.path.join(d, f'{name}.pdb')); t.save_dcd(os.path.join(d, f'{name}.dcd'))
        done.append(name); print(f'  superposed {name}: {t.n_frames} frames, core CA RMSD of frame 0 to the reference {rms:.2f} A', flush=True)
    return done


def pymol_run(spec, tag):
    p = os.path.join(W_, f'spec_{tag}.json'); json.dump(spec, open(p, 'w'))
    r = subprocess.run([PYMOL, '-cq', SCENE, '--', p], capture_output=True, text=True)
    if r.returncode != 0 or 'Error' in r.stdout + r.stderr:
        print(r.stdout[-2000:], r.stderr[-2000:])
    assert r.returncode == 0


def panel(name, style=None):
    sec, _, A, B, iso, sty = CL[name]; d = os.path.join(a.out, sec, 'pymol')
    p = dict(name=name, pdb=os.path.join(d, f'{name}.pdb'), dcd=os.path.join(d, f'{name}.dcd'), style=style or sty, alpha=A, partner=B, iso=iso)
    if p['style'] == 'unfold':
        p['units'] = UNITS
    return p


AROM = '(resi 449 and resn PHE or resi 451 and resn TYR) and name CG+CD1+CD2+CE1+CE2+CZ'
VIEWS = {'loop': dict(view_sel='chain A and resi 30-47', buffer=9, look=['chain A and name CA', 'chain A and resi 32-45'],
                      up=['chain A and resi 1-200 and name CA', 'chain A and resi 260-400 and name CA']),
         'groove': dict(view_sel='(chain B and resi 187+192+203+223) or (chain A and resi 444-451)', buffer=8,
                        look=['chain B and resi 160-244 and name CA', f'chain A and {AROM}'], up=['chain A and name CA', 'chain B and resi 160-244 and name CA']),
         'pocket': dict(view_sel='chain A and resi 108-120', buffer=5, look=['chain A and name CA', 'chain A and resi 110-118 and name CA'],
                        up=['chain A and resi 1-200 and name CA', 'chain A and resi 260-400 and name CA']),
         'unfold': dict(view_sel='chain A and name CA', buffer=1),
         'overview': dict(view_sel='chain A and name CA', buffer=2)}


def view_for(style, alpha):
    v = json.loads(json.dumps(VIEWS[style]))
    if alpha != 'A':                                   # 9M1L: alpha is chain B, TBCE chain A (loop views only name chain A)
        assert style == 'loop', style
        fix = lambda x: x.replace('chain A', f'chain {alpha}')
        v = {k: (fix(x) if isinstance(x, str) else [fix(y) for y in x] if isinstance(x, list) else x) for k, x in v.items()}
    return v


def session(name):
    sec, _, A, B, iso, sty = CL[name]
    spec = dict(panels=[panel(name)], mode='pse', out=os.path.join(a.out, sec, 'pymol', f'{name}.pse'), **view_for(sty, A))
    pymol_run(spec, f'pse_{name}')


# ------------------------------------------------------------------ read-out strips under each movie panel
def strip(name, kind, width_px, height_px, path):
    r = series(name); dt = r['ns_per_frame']; dpi = 100
    fig = plt.figure(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)
    ax = fig.add_axes([0.08, 0.30, 0.88, 0.62])
    c = ISOC[CL[name][4]]
    w = max(1, int(round(1.0 / dt)))
    rmean = lambda x: np.convolve(x, np.ones(w), 'same') / np.convolve(np.ones(len(x)), np.ones(w), 'same')
    if kind == 'loop':
        q = np.array(r['loop32_45']['Q_lam1.2']); rm = np.array(r['loop32_45']['rmsd_loop_A']); t = np.arange(len(q)) * dt
        st = np.where((q >= QN) & (rm <= RN), 0, np.where(q < QU, 2, 1))
        ax.imshow(st[None, :], extent=[0, t[-1], -0.02, 0.08], aspect='auto', cmap=mpl.colors.ListedColormap(['#2b2b2a', '#b9b8b0', '#f2f1ec']), vmin=0, vmax=2)
        ax.plot(t, rmean(q), color=c, lw=1.6); ax.axhline(QN, color=INK2, lw=0.8, ls='--')
        ax.set_ylim(-0.02, 0.85); ax.set_ylabel('loop Q', fontsize=9); ax.text(t[-1], 0.1, 'state: N dark / I grey / U light', ha='right', fontsize=7.5, color=INK2)
    elif kind == 'groove':
        gg = r['groove']; dmin = np.minimum(np.array(gg['d_ring_PHE187_A']), np.array(gg['d_ring_TRP192_A'])); t = np.arange(len(dmin)) * dt
        ax.plot(t, rmean(dmin), color=c, lw=1.6); ax.axhline(6.5, color=INK2, lw=0.8, ls='--')
        ax.set_ylim(3, 14); ax.set_ylabel(f'{gg["tail_residue"]} ring to\nF187/W192 (Å)', fontsize=8)
    elif kind == 'unfold':
        G = ['loop32_45', 'unit1', 'unit2', 'unit3', 'unit0', 'seg110_120', 'unit5', 'unit4']
        M = np.array([rmean(np.array(r['units'][g])) for g in G]); t = np.arange(M.shape[1]) * dt
        ax.imshow(M, aspect='auto', extent=[0, t[-1], len(G) - 0.5, -0.5], cmap='Blues', vmin=0.2, vmax=0.9)
        ax.set_yticks(range(len(G))); ax.set_yticklabels(['loop', 'u1', 'u2', 'u3', 'u0', '110–120', 'u5', 'u4'], fontsize=6.5)
    elif kind == 'pocket':
        res = r['res114']
        if 'd_OG114_O110_A' in res:
            d = np.array(res['d_OG114_O110_A']); t = np.arange(len(d)) * dt
            ax.plot(t, rmean(d), color=c, lw=1.6); ax.axhline(3.5, color=INK2, lw=0.8, ls='--'); ax.set_ylim(2, 8)
            ax.set_ylabel('Ser114 OG –\nVal110 O (Å)', fontsize=8)
        else:
            n = np.array(res['n_contacts_sc_4.5A'], float); t = np.arange(len(n)) * dt
            ax.plot(t, rmean(n), color=c, lw=1.6); ax.set_ylabel('Ile114 side-chain\ncontacts', fontsize=8)
    ax.set_xlim(0, t[-1]); ax.tick_params(labelsize=8); ax.set_xlabel('time (ns)', fontsize=8, labelpad=1)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    fig.savefig(path, dpi=dpi, facecolor='white'); plt.close(fig)
    x0 = ax.get_position().x0 * width_px; x1 = ax.get_position().x1 * width_px
    y0 = (1 - ax.get_position().y1) * height_px; y1 = (1 - ax.get_position().y0) * height_px
    return dict(x0=x0, x1=x1, y0=y0, y1=y1, tmax=float(t[-1]))


def movie(tag, names, kind, title, sub, style=None, width=1280, height=760):
    sec = CL[names[0]][0]; out = os.path.join(a.out, sec, f'movie_{tag}.mp4')
    fr = os.path.join(W_, f'frames_{tag}'); shutil.rmtree(fr, ignore_errors=True); os.makedirs(fr)
    sty = style or CL[names[0]][5]
    spec = dict(panels=[panel(n, sty) for n in names], mode='frames', out=fr, width=width, height=height, smooth=5, **view_for(sty, CL[names[0]][2]))
    pymol_run(spec, f'mov_{tag}')
    SH, TH_ = 150, 92
    strips = [strip(n, kind, width, SH, os.path.join(fr, f'strip_{n}.png')) for n in names]
    times = [series(n)['clip']['clip_time_ns'] for n in names]
    nfr = min(len(x) for x in times)
    fB = ImageFont.truetype(FONTB, 30); fR = ImageFont.truetype(FONT, 21); fP = ImageFont.truetype(FONTB, 24); fC = ImageFont.truetype(FONT, 24)
    comp = os.path.join(fr, 'comp'); os.makedirs(comp)
    labels = [n.replace('_', ' ') for n in names]
    for k in range(nfr):
        img = Image.new('RGB', (width * len(names), TH_ + height + SH), 'white'); dr = ImageDraw.Draw(img)
        dr.text((24, 14), title, font=fB, fill=INK); dr.text((24, 54), sub, font=fR, fill=INK2)
        for j, n in enumerate(names):
            x = j * width
            im = Image.open(os.path.join(fr, n, f'f{k + 1:04d}.png')).convert('RGB'); img.paste(im, (x, TH_))
            dr.text((x + 20, TH_ + 10), labels[j], font=fP, fill=ISOC[CL[n][4]])
            dr.text((x + width - 190, TH_ + 10), f't = {times[j][k]:6.1f} ns', font=fC, fill=INK)
            s = Image.open(os.path.join(fr, f'strip_{n}.png')).convert('RGB'); img.paste(s, (x, TH_ + height))
            g = strips[j]; cx = g['x0'] + (g['x1'] - g['x0']) * times[j][k] / g['tmax']
            dr.line([(x + cx, TH_ + height + g['y0']), (x + cx, TH_ + height + g['y1'])], fill='#e34948', width=3)
            if j:
                dr.line([(x, TH_), (x, TH_ + height + SH)], fill=GRID, width=2)
        img.save(os.path.join(comp, f'c{k:04d}.png'))
    subprocess.run([FFMPEG, '-y', '-loglevel', 'error', '-framerate', '25', '-i', os.path.join(comp, 'c%04d.png'),
                    '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', '-crf', '18', '-movflags', '+faststart', out], check=True)
    print('  movie', out, nfr, 'frames', flush=True)


def have(*names):
    return all(os.path.exists(os.path.join(a.out, CL[n][0], 'pymol', f'{n}.dcd')) for n in names)


print('superposing clips on the common reference')
if not a.skip_prep:
    done = superpose_all()
    for n in done:
        session(n); print('  session', n, flush=True)
    if a.prep_names:
        sys.exit(0)
MOV = [('free_TUBA1A_r1_loop', ['free_TUBA1A_r1'], 'loop', 'Free TUBA1A monomer (r1, 200 ns, 310 K): loop 32–45 leaves its starting conformation',
        'blue: loop 32–45, sticks: residues 39/42/43/45; bottom: loop Q (1-ns mean) and state; cartoon lightly smoothed for display'),
       ('heterodimer_TUBA1A_vs_TUBA8_loop', ['heterodimer_TUBA1A_r1', 'heterodimer_TUBA8_r0'], 'loop',
        'αβ heterodimer (9M1M), 50 ns: TUBA1A loop stays near its start, TUBA8 loop does not', 'β-tubulin in light blue-grey; same camera for both; bottom: loop Q and N/I/U state'),
       ('TBCE9M1L_TUBA1A_vs_TUBA8_loop', ['TBCE9M1L_TUBA1A_r0', 'TBCE9M1L_TUBA8_r1'], 'loop',
        'TBCE complex (9M1L), 50 ns: the TUBA8 loop opens up (U 95 % of r1)', 'TBCE in light blue-grey; bottom: loop Q and N/I/U state'),
       ('TBCB_TUBA8_r1_groove', ['TBCB_TUBA8_r1'], 'groove', 'TBCB CAP-Gly groove holds the TUBA8 tail end (Phe449), 50 ns',
        'grey sticks: TBCB F187, K188, W192, K203, F223; orange: α tail 438–449; bottom: ring–ring distance (dashed 6.5 Å)'),
       ('TBCB_TUBA1A_r2_groove', ['TBCB_TUBA1A_r2'], 'groove', 'TBCB CAP-Gly groove holds the TUBA1A tail end (Tyr451), 50 ns',
        'grey sticks: TBCB F187, K188, W192, K203, F223; blue: α tail 438–451; bottom: ring–ring distance (dashed 6.5 Å)'),
       ('450K_TUBA1A_vs_TUBA8_unfolding', ['md450_TUBA1A_r0', 'md450_TUBA8_r0'], 'unfold', '450 K unfolding, 20 ns (r0 of each): the 32–45 loop goes first',
        'colours: folding units 0–5 (blue, orange, aqua, yellow, pink, green), loop 32–45 violet, 110–120 red; bottom: Q per unit'),
       ('450K_TUBA8_r0_unfolding', ['md450_TUBA8_r0'], 'unfold', 'TUBA8 at 450 K (r0, 20 ns): the 32–45 loop goes first, units 4 and 5 last',
        'colours: folding units 0–5 (blue, orange, aqua, yellow, pink, green), loop 32–45 violet, 110–120 red; bottom: Q per unit (1-ns mean)'),
       ('Ser114_vs_Ile114_pocket', ['free_TUBA8_r0', 'free_TUBA1A_r0'], 'pocket', 'Residue 114 pocket: Ser114 (TUBA8) vs Ile114 (TUBA1A), free monomer, 200 ns',
        'sticks: 114 and every residue within 4.5 Å; dashed: Ser114 OG – Val110 O; bottom: that distance (TUBA8) / side-chain contacts (TUBA1A)')]
only = set(x for x in a.only.split(',') if x)
if not a.no_movies and not a.stills_only:
    for tag, names, kind, title, sub in MOV:
        if only and tag not in only:
            continue
        if not have(*names):
            print('  skip movie (clip missing):', tag); continue
        movie(tag, names, kind, title, sub, style='pocket' if kind == 'pocket' else None)


# ------------------------------------------------------------------ stills: representative frames chosen by a stated rule
def at_clip(name, key_fn):
    r = series(name); idx = np.array(r['clip']['clip_frame_index']); tt = np.array(r['clip']['clip_time_ns'])
    return r, idx, tt, np.asarray(key_fn(r), float)[idx]


def closest(values, target, mask):
    i = np.where(mask)[0]; return int(i[np.argmin(np.abs(values[i] - target))])


def compose(tiles, captions, title, sub, path, ncol=None):
    ims = [Image.open(t).convert('RGB') for t in tiles]; w, h = ims[0].size; ncol = ncol or len(ims); nrow = int(np.ceil(len(ims) / ncol))
    fB = ImageFont.truetype(FONTB, 40); fR = ImageFont.truetype(FONT, 28); fC = ImageFont.truetype(FONTB, 32)
    canvas = Image.new('RGB', (w * ncol, 120 + nrow * (h + 60)), 'white'); dr = ImageDraw.Draw(canvas)
    dr.text((30, 18), title, font=fB, fill=INK); dr.text((30, 70), sub, font=fR, fill=INK2)
    for k, (im, (cap, col)) in enumerate(zip(ims, captions)):
        x, y = (k % ncol) * w, 120 + (k // ncol) * (h + 60)
        dr.text((x + 24, y + 8), cap, font=fC, fill=col); canvas.paste(im, (x, y + 56))
    canvas.save(path, dpi=(300, 300)); print('  still', path, flush=True)


def still(names, states, style, tag, labels=None, width=1500, height=1100):
    d = os.path.join(W_, f'still_{tag}'); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    spec = dict(panels=[panel(n, style) for n in names], mode='still', states=states, out=d, width=width, height=height, labels=labels,
                per_panel_view=True, **view_for(style, CL[names[0]][2]))
    pymol_run(spec, f'still_{tag}')
    return [[os.path.join(d, f'{n}_s{s_:04d}.png') for s_ in st] for n, st in zip(names, states)]


def make_stills():
    after10 = lambda tt: tt >= 10.0
    if have('free_TUBA1A_r1'):
        r, idx, tt, q = at_clip('free_TUBA1A_r1', lambda r: r['loop32_45']['Q_lam1.2'])
        rm = np.asarray(r['loop32_45']['rmsd_loop_A'], float)[idx]
        st = np.where((q >= QN) & (rm <= RN), 0, np.where(q < QU, 2, 1))
        pick = []
        for s_ in (0, 1, 2):
            m = st == s_
            if m.any(): pick.append((s_, closest(q, np.median(q[m]), m)))
        tiles = still(['free_TUBA1A_r1'], [[k + 1 for _, k in pick]], 'loop', 'loop_states')[0]
        caps = [(f'{"NIU"[s_]}: t = {tt[k]:.1f} ns, Q = {q[k]:.2f}, loop RMSD = {rm[k]:.1f} Å', BLUE) for s_, k in pick]
        compose(tiles, caps, 'Loop 32–45 states in one free TUBA1A trajectory (r1)', 'each panel: the clip frame whose Q is closest to the median Q of that state; same camera',
                os.path.join(a.out, '04_loop32-45_states', 'still_loop_N_I_U_TUBA1A_r1.png'))
    if have('free_TUBA8_r0', 'free_TUBA1A_r0'):
        r8, i8, t8, d8 = at_clip('free_TUBA8_r0', lambda r: r['res114']['d_OG114_O110_A'])
        r1, i1, t1, n1 = at_clip('free_TUBA1A_r0', lambda r: r['res114']['n_contacts_sc_4.5A'])
        k8 = closest(d8, np.median(d8[after10(t8)]), after10(t8)); k1 = closest(n1, np.median(n1[after10(t1)]), after10(t1))
        tl = still(['free_TUBA8_r0', 'free_TUBA1A_r0'], [[k8 + 1], [k1 + 1]], 'pocket', 'pocket',
                   labels=[('resi 114 and name CB', '114'), ('resi 110 and name CA', 'V110')])
        compose([tl[0][0], tl[1][0]], [(f'TUBA8 Ser114: t = {t8[k8]:.1f} ns, OG–O(V110) = {d8[k8]:.1f} Å (median frame)', ORANGE),
                                       (f'TUBA1A Ile114: t = {t1[k1]:.1f} ns, {n1[k1]:.0f} side-chain contacts (median frame)', BLUE)],
                'Residue 114 pocket, free monomers', 'sticks: residue 114 and all residues within 4.5 Å of its side chain; frames at the median of the quantity shown',
                os.path.join(a.out, '01_Ser114_folding_cost', 'still_Ser114_vs_Ile114_pocket.png'))
    if have('TBCB_TUBA8_r1', 'TBCB_TUBA1A_r2'):
        tiles, caps = [], []
        for n, col in (('TBCB_TUBA8_r1', ORANGE), ('TBCB_TUBA1A_r2', BLUE)):
            r, idx, tt, dmin = at_clip(n, lambda r: np.minimum(np.array(r['groove']['d_ring_PHE187_A']), np.array(r['groove']['d_ring_TRP192_A'])))
            k = closest(dmin, np.median(dmin[after10(tt)]), after10(tt))
            tl = still([n], [[k + 1]], 'groove', f'groove_{n}', labels=[('chain B and resi 187 and name CZ', 'F187'), ('chain B and resi 192 and name CZ2', 'W192'),
                                                                     ('chain B and resi 223 and name CZ', 'F223'), ('chain B and resi 203 and name NZ', 'K203'),
                                                                     ('chain B and resi 188 and name NZ', 'K188'), (f'chain A and {AROM.split(" and name")[0]} and name CZ', series(n)['groove']['tail_residue'])])
            tiles.append(tl[0][0]); caps.append((f'{n.replace("_", " ")}: t = {tt[k]:.1f} ns, ring distance {dmin[k]:.1f} Å (median frame)', col))
        compose(tiles, caps, 'TBCB CAP-Gly groove around the tail-end aromatic', 'grey sticks: TBCB groove residues; coloured: α tail; frames at the median ring–ring distance after 10 ns',
                os.path.join(a.out, '02_TBCB_tail_clamp', 'still_TBCB_groove_TUBA8_vs_TUBA1A.png'))
    for n, col in (('md450_TUBA8_r0', ORANGE), ('md450_TUBA1A_r0', BLUE)):
        if not have(n): continue
        tt = np.array(series(n)['clip']['clip_time_ns']); ks = [int(np.argmin(np.abs(tt - x))) for x in (0, 2, 5, 10, 20)]
        tl = still([n], [[k + 1 for k in ks]], 'unfold', f'unfold_{n}', width=1100, height=1000)[0]
        compose(tl, [(f't = {tt[k]:.1f} ns', col) for k in ks], f'{n.replace("md450_", "").replace("_r0", "")} at 450 K (r0): snapshots',
                'colours: units 0–5 blue, orange, aqua, yellow, pink, green; loop 32–45 violet; 110–120 red', os.path.join(a.out, '06_450K_unfolding_order', f'still_450K_{n}_snapshots.png'))
    ov = 'free_TUBA8_r0' if have('free_TUBA8_r0') else 'free_TUBA1A_r1'
    tl = still([ov], [[1]], 'overview', 'overview', labels=[('resi 39 and name CA', 'loop 32-45'), ('resi 114 and name CB', '114'),
                                                            ('resn GTP and name PG', 'GTP'), ('resi 437 and name CA', 'C-terminus 437')], width=1800, height=1300)[0]
    compose(tl, [(f'{ov.split("_")[1]} free monomer, frame 0', ORANGE if 'TUBA8' in ov else BLUE)], 'α-tubulin: the sites in this deck',
            'loop 32–45 with side chains 39/42/43/45; residue 114 as spheres (Ser in TUBA8, Ile elsewhere); GTP + Mg2+.\n'
            'The free-monomer builds end at residue 437; the tail 438–449/451 that TBCB and TBCE bind is modelled only in the complexes.',
            os.path.join(a.out, '00_overview', 'still_alpha_overview_sites.png'))


if not a.no_movies or a.stills_only:
    pass
if a.stills_only or not a.only or a.with_stills:
    os.makedirs(os.path.join(a.out, '00_overview'), exist_ok=True)
    make_stills()
