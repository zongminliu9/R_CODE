"""Lab-meeting figures (PNG 300 dpi + vector PDF) from the compiled tables and the a31 clip series.

Every number is read from the outputs of compile_interim_tables.py (interim/*/tables_*.json), the loop CV files
(interim/*/loopcv_md310_*.json), the PMF / tauRAMD result files and the a31 series; values the reports quote are
re-derived here and asserted, so a figure cannot drift from the report silently.
Usage: a32_labmeeting_figures.py --out <labmeeting dir> --work <labmeeting_work dir>
"""
import os, json, glob, argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser(); ap.add_argument('--out', required=True); ap.add_argument('--work', required=True)
a = ap.parse_args()
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
I = lambda *p: os.path.join(ROOT, 'interim', *p)
J = lambda p: json.load(open(p))

BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
GRAY, INK, INK2, GRID, AXIS = '#898781', '#0b0b0b', '#52514e', '#e1e0d9', '#c3c2b7'
ISO = {'TUBA1A': BLUE, 'TUBA8': ORANGE}
mpl.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 12,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold',
    'axes.titlesize': 13, 'axes.titleweight': 'bold', 'axes.labelsize': 12, 'axes.edgecolor': AXIS,
    'axes.labelcolor': INK, 'axes.titlecolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
    'grid.linewidth': 0.6, 'axes.axisbelow': True, 'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'savefig.dpi': 300, 'legend.frameon': False, 'lines.linewidth': 2, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
KT = 0.616
made = []


def save(fig, sub, name):
    d = os.path.join(a.out, sub); os.makedirs(d, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(d, f'{name}.{ext}'), bbox_inches='tight', facecolor='white')
    plt.close(fig); made.append(f'{sub}/{name}')


def rmean(x, w):
    # centred running mean whose edge windows are averaged over the points they actually contain (no zero padding)
    x = np.asarray(x, float); k = np.ones(w)
    return np.convolve(x, k, mode='same') / np.convolve(np.ones_like(x), k, mode='same')


def fv(v, spec='+.2f'):
    return format(v, spec).replace('-', '\u2212')


def close(x, y, tol, what):
    assert abs(x - y) <= tol, f'{what}: {x} vs report {y}'


T6, T7, T8, TM = (J(I('20261106_1200', 'tables_20261106_1200.json')), J(I('20261107_1730', 'tables_20261107_1730.json')),
                  J(I('20261108_final', 'tables_20261108_final.json')), J(I('20261110_m', 'tables_20261110_m.json')))

# ---------------------------------------------------------------- 1. folding free energies
dd8, ddm = T8['ddG']['this_pass'], TM['ddG']['this_pass']
fepm = TM['fep']['this_pass']
i114 = dd8['a1a_114SER  [ref a1a_114SER_unfolded]']; s114 = dd8['a8_114ILE  [ref a8_114ILE_unfolded]']
# loop-swapped backgrounds: folded leg of the swapped build minus the native unfolded pentapeptide leg
sw1 = fepm['a1a_loop8_114SER_folded']['summary']; sw8 = fepm['a8_loop1A_114ILE_folded']['summary']
sw1_dd = sw1['mean_dG_kcal'] - i114['unfolded']['mean_dG_kcal']; sw1_sd = np.hypot(sw1['sd_dG_kcal'], i114['unfolded']['sd_dG_kcal'])
sw8_dd = sw8['mean_dG_kcal'] - s114['unfolded']['mean_dG_kcal']; sw8_sd = np.hypot(sw8['sd_dG_kcal'], s114['unfolded']['sd_dG_kcal'])
close(sw1_dd, 4.50, 0.006, 'loop-swap I114S'); close(sw8_dd, -4.34, 0.006, 'loop-swap S114I')
close(i114['ddG_fold_kcal'], 4.04, 0.006, 'I114S'); close(s114['ddG_fold_kcal'], -4.10, 0.006, 'S114I')
rows = [('I114S in TUBA1A', i114['ddG_fold_kcal'], i114['sd_kcal'], BLUE, True, '3 + 3 replicas'),
        ('S114I in TUBA8', s114['ddG_fold_kcal'], s114['sd_kcal'], ORANGE, True, '3 + 3; folded legs: D-Gln39/D-Asp45*'),
        ('I114S in TUBA1A, TUBA8 loop', sw1_dd, sw1_sd, BLUE, True, 'folded 2/3 done'),
        ('S114I in TUBA8, TUBA1A loop', sw8_dd, sw8_sd, ORANGE, True, 'folded 3/3'),
        ('I110V in TUBA1A', ddm['a1a_110VAL  [ref a1a_110VAL_unfolded]']['ddG_fold_kcal'], ddm['a1a_110VAL  [ref a1a_110VAL_unfolded]']['sd_kcal'], BLUE, True, '2 + 2'),
        ('A294S in TUBA1A', ddm['a1a_294SER  [ref a1a_294SER_unfolded]']['ddG_fold_kcal'], ddm['a1a_294SER  [ref a1a_294SER_unfolded]']['sd_kcal'], BLUE, True, '2 + 2'),
        ('A270V in TUBA1A', ddm['a1a_270VAL  [ref a1a_270VAL_unfolded]']['ddG_fold_kcal'], ddm['a1a_270VAL  [ref a1a_270VAL_unfolded]']['sd_kcal'], BLUE, False, 'not converged')]
close(rows[4][1], 0.10, 0.006, 'I110V'); close(rows[5][1], 0.19, 0.006, 'A294S'); close(rows[6][1], -1.49, 0.006, 'A270V')
fig, ax = plt.subplots(figsize=(9.5, 4.8))
y = np.arange(len(rows))[::-1]
for yi, (lab, v, sd, c, filled, note) in zip(y, rows):
    ax.errorbar(v, yi, xerr=sd, fmt='o', ms=10, color=c, mfc=c if filled else 'white', mew=2, elinewidth=2, capsize=0, zorder=3)
    ax.text(v, yi + 0.3, f'{fv(v)} ± {sd:.2f}', va='bottom', ha='center', color=INK, fontsize=10.5)
    ax.text(1.01, yi, note, va='center', ha='left', color=INK2, fontsize=10, transform=ax.get_yaxis_transform())
ax.axvline(0, color=AXIS, lw=1, zorder=1)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.set_xlim(-6.2, 6.2); ax.set_ylim(-0.6, len(rows) - 0.3); ax.grid(axis='y', visible=False)
ax.set_xlabel('ΔΔG$_{fold}$ (kcal/mol)   > 0: mutation destabilises the folded α')
ax.set_title('Only the Ser/Ile switch at 114 costs folding stability (≈ 4 kcal/mol), in either loop background', loc='left')
ax.text(0.0, -0.17, '* folded legs started before the chirality fix: Gln39 and Asp45 are D (29–34 Å from 114); the all-L runs (I114S in TUBA1A, S114I with the TUBA1A loop) agree',
        transform=ax.transAxes, fontsize=9, color=INK2)
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, ms=9, label='TUBA1A background'),
                   Line2D([], [], marker='o', ls='', color=ORANGE, ms=9, label='TUBA8 background'),
                   Line2D([], [], marker='o', ls='', color=GRAY, mfc='white', mew=2, ms=9, label='open = not converged')],
          loc='lower left', fontsize=10)
save(fig, '01_Ser114_folding_cost', 'fig1a_ddG_fold_all_sites')

# alchemical path: cumulative dG along the five physical steps, folded vs unfolded leg, per replica
a27 = T8['a27']['this_pass']; ORDER = ['new torsions on', 'old charges off', 'sterics swap', 'new charges on', 'old torsions off']
pdec = T8['path_decomposition']['this_pass']
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={'width_ratios': [1, 1, 1.05]})
for ax, (mut, base, c, title) in zip(axs[:2], [('I114S', 'a1a_114SER', BLUE, 'I114S in TUBA1A (Ile → Ser)'), ('S114I', 'a8_114ILE', ORANGE, 'S114I in TUBA8 (Ser → Ile)')]):
    for leg, ls in (('folded', '-'), ('unfolded', '--')):
        curves = []
        for rep, rr in a27[f'{base}_{leg}'].items():
            cum = np.concatenate([[0], np.cumsum([rr['steps'][s] for s in ORDER])]); curves.append(cum)
            ax.plot(range(6), cum, ls, color=c, lw=1, alpha=0.35)
        m = np.mean(curves, axis=0); ax.plot(range(6), m, ls, color=c, lw=2.5, marker='o', ms=6, label=f'{leg} leg (mean, thin = replicas)')
        ax.text(5.1, m[-1], f'{leg}\n{fv(m[-1], "+.1f")}', va='center', fontsize=10, color=INK)
    ax.axhline(0, color=AXIS, lw=1)
    ax.set_xticks(range(6)); ax.set_xticklabels(['start'] + [s.replace(' ', '\n') for s in ORDER], fontsize=9)
    ax.set_ylabel('cumulative ΔG (kcal/mol)'); ax.set_title(title, loc='left'); ax.set_xlim(-0.2, 5.9)
    ax.legend(loc='lower left' if mut == 'I114S' else 'lower right', fontsize=9)
ax = axs[2]; w = 0.38
# group the steps by what they physically switch, in the Ile -> Ser direction for both calculations:
# I114S runs Ile -> Ser (old = Ile, new = Ser); S114I runs Ser -> Ile, so its steps are mapped and sign-reversed
PHYS = [('Ser OH\ncharges\non', 'new charges on', 'old charges off'), ('packing\n(sterics\nswap)', 'sterics swap', 'sterics swap'),
        ('Ile\ncharges\noff', 'old charges off', 'new charges on'), ('Ile\ntorsions\noff', 'old torsions off', 'new torsions on'),
        ('Ser\ntorsions\non', 'new torsions on', 'old torsions off')]
x = np.arange(len(PHYS))
dI = [pdec['I114S'][a_] for _, a_, _b in PHYS]; dS = [-pdec['S114I（暂定）'][b_] for _, _a, b_ in PHYS]
close(sum(dI), i114['ddG_fold_kcal'], 0.02, 'path decomposition sum')
ax.bar(x - w / 2 - 0.01, dI, w, color=BLUE, label=f'I114S in TUBA1A (sum {fv(sum(dI))})')
ax.bar(x + w / 2 + 0.01, dS, w, color=ORANGE, label=f'−(S114I in TUBA8) (sum {fv(sum(dS))}, provisional)')
for xs_, vals in ((x - w / 2, dI), (x + w / 2, dS)):
    for xi, v in zip(xs_, vals):
        if abs(v) > 0.3: ax.text(xi, v + (0.06 if v > 0 else -0.06), fv(v, '+.1f'), ha='center', va='bottom' if v > 0 else 'top', fontsize=9, color=INK)
ax.axhline(0, color=AXIS, lw=1); ax.set_xticks(x); ax.set_xticklabels([q_[0] for q_ in PHYS], fontsize=9)
ax.set_ylabel('share of ΔΔG(Ile → Ser)  (kcal/mol)'); ax.set_title('Where the 4 kcal/mol comes from', loc='left'); ax.grid(axis='x', visible=False)
ax.set_ylim(-0.6, 3.0); ax.legend(loc='upper right', fontsize=9)
fig.suptitle('Alchemical path at residue 114: packing (sterics swap) + buried hydroxyl (charges on)', x=0.01, ha='left', fontsize=14, fontweight='bold', y=1.03)
save(fig, '01_Ser114_folding_cost', 'fig1b_alchemical_path_114')

# convergence: reverse-cumulative ddG (last x of the data) and replica blocks
rc = T8['ddG_reverse_cumulative']['this_pass']
fig, axs = plt.subplots(1, 2, figsize=(13, 4.4))
ax = axs[0]
for key, c, lab in (('I114S', BLUE, 'I114S'), ('S114I（暂定）', ORANGE, 'S114I')):
    xs = sorted(rc[key], key=float); ax.plot([float(k) * 100 for k in xs], [rc[key][k] for k in xs], '-o', color=c, ms=7, label=lab)
    ax.text(102, rc[key]['1.0'], f'{lab} {fv(rc[key]["1.0"])}', va='center', fontsize=10, color=INK)
ax.axhline(0, color=AXIS, lw=1); ax.set_xlim(15, 125); ax.set_ylim(-5.5, 5.5)
ax.set_xlabel('data used, counted back from the end (%)'); ax.set_ylabel('ΔΔG$_{fold}$ (kcal/mol)')
ax.set_title('Estimate is flat in the amount of data used', loc='left'); ax.legend(loc='center left', fontsize=10)
ax = axs[1]
pos = 0; ticks = []
for base, c in (('a1a_114SER', BLUE), ('a8_114ILE', ORANGE)):
    for leg in ('folded', 'unfolded'):
        for rep, rr in a27[f'{base}_{leg}'].items():
            b = rr['blocks']; ax.plot(pos + np.linspace(-0.3, 0.3, len(b)), b, '-o', color=c, ms=5, lw=1.5)
            ticks.append((pos, rep)); pos += 1
        ax.text(pos - 2, -0.13, f'{"I114S" if base.startswith("a1a") else "S114I"} {leg}', ha='center', va='top', fontsize=9.5, color=INK, transform=ax.get_xaxis_transform())
        pos += 0.6
ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks], fontsize=9)
ax.set_ylabel('ΔG of the leg, 4 consecutive blocks'); ax.set_title('Within-replica drift (4 blocks per replica)', loc='left'); ax.grid(axis='x', visible=False)
save(fig, '01_Ser114_folding_cost', 'fig1c_fep_convergence_114')

# ---------------------------------------------------------------- 2. tail binding: ddG(Tyr -> Phe), tauRAMD
tr = T8['tail_rbfe']
def g(k): return tr[k]['value'], tr[k]['sd']
rowsB = [('TBCB', 'TUBA1A tail (Y451F)', *g('TBCB_Y451F_Tuba1a_tail'), BLUE, True),
         ('TBCB', 'TUBA8 tail (−F449Y)', -g('TBCB_F449Y_Tuba8_tail')[0], g('TBCB_F449Y_Tuba8_tail')[1], ORANGE, True),
         ('TBCB-K188A', 'TUBA1A tail (Y451F)', *g('TBCB-K188A_Y451F_Tuba1a_tail'), BLUE, True),
         ('TBCB-K188A', 'TUBA8 tail (−F449Y)', -g('TBCB-K188A_F449Y_Tuba8_tail')[0], g('TBCB-K188A_F449Y_Tuba8_tail')[1], ORANGE, True),
         ('TBCE', 'TUBA1A tail (Y451F)', *g('TBCE_Y451F_Tuba1a_tail'), BLUE, False),
         ('TBCE', 'TUBA8 tail (−F449Y)', -g('TBCE_F449Y_Tuba8_tail')[0], g('TBCE_F449Y_Tuba8_tail')[1], ORANGE, True)]
close(rowsB[0][2], 0.72, 0.006, 'TBCB Y451F'); close(rowsB[1][2], 0.74, 0.006, 'TBCB F449Y (sign flipped)')
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={'width_ratios': [1.15, 1]})
ax = axs[0]; ys = []
yy = 0
for i, (rec, lab, v, sd, c, ok) in enumerate(rowsB):
    if i and rowsB[i - 1][0] != rec: yy += 0.6
    ax.errorbar(v, -yy, xerr=sd, fmt='o', ms=10, color=c, mfc=c if ok else 'white', mew=2, elinewidth=2)
    ax.text(v + sd + 0.08, -yy, f'{fv(v)} ± {sd:.2f}' + ('' if ok else '  unresolved'), va='center', fontsize=10, color=INK)
    ys.append((-yy, f'{rec}: {lab}')); yy += 1
ax.axvline(0, color=AXIS, lw=1); ax.set_yticks([q[0] for q in ys]); ax.set_yticklabels([q[1] for q in ys], fontsize=10)
ax.set_xlim(-1.4, 2.0); ax.grid(axis='y', visible=False)
ax.set_xlabel('ΔΔG$_{bind}$(Tyr → Phe at the tail end), kcal/mol   > 0: Tyr binds better')
cI, cE = tr['coupling_K188_hydroxyl_Tuba1a_tail'], tr['coupling_K188_hydroxyl_Tuba8_tail']
ax.set_title('TBCB prefers a Tyr-ended α tail by ≈ 0.7 kcal/mol', loc='left')
ax.text(0.02, 0.56, f'K188 × hydroxyl coupling:\n{fv(cI["value"])} ± {cI["sd"]:.2f} (TUBA1A tail)\n{fv(cE["value"])} ± {cE["sd"]:.2f} (TUBA8 tail)',
        fontsize=9.5, color=INK2, transform=ax.transAxes, va='top')
# tauRAMD survival curves from the raw dissociation times
RAMD = os.path.join(a.work, 'data', 'runs', 'ramd'); ramd_tab = T8['ramd']
def surv(sid):
    t = [];
    for f in sorted(glob.glob(os.path.join(RAMD, sid, 's*', 'ramd_times.json'))):
        # a run that never dissociated is censored at its length, 1e7 steps x 2 fs = 20 ns (as a22 counts it)
        t += [(20.0 if r['t_ns'] is None else r['t_ns'], r['censored'] or r['t_ns'] is None) for r in J(f)]
    t.sort(); n = len(t); S = [1.0]; tt = [0.0]; alive = n
    for ti, cen in t:
        if not cen: tt += [ti, ti]; S += [S[-1], S[-1] * (alive - 1) / alive]
        alive -= 1
    return np.array(tt), np.array(S), n, float(np.median([x[0] for x in t]))
ax = axs[1]
for sid, c, ls, lab in (('t_a8_cplx', ORANGE, '--', 'TUBA8 tail, Phe449 (native)'), ('t_a8FY_cplx', ORANGE, '-', 'TUBA8 tail, F449Y'),
                        ('t_a8_K188A_cplx', GRAY, '--', 'K188A, Phe449'), ('t_a8FY_K188A_cplx', GRAY, '-', 'K188A, F449Y')):
    tt, S, n, med = surv(sid); assert n == ramd_tab[sid]['n'] == 30
    close(med, ramd_tab[sid]['tau_ns'], 1e-6, f'tauRAMD median {sid}')
    ax.step(np.r_[tt, 12], np.r_[S, S[-1]], where='post', color=c, ls=ls, lw=2.2, label=f'{lab}: τ = {med:.2f} ns')
ax.axhline(0.5, color=AXIS, lw=1, ls=':'); ax.set_xlim(0, 10); ax.set_ylim(0, 1.02)
ax.set_xlabel('time under random force (ns)'); ax.set_ylabel('fraction of 30 trajectories still bound')
ax.set_title('τRAMD, TBCB + TUBA8 tail (30 runs each)', loc='left'); ax.legend(fontsize=9.5, loc='upper right')
save(fig, '02_TBCB_tail_clamp', 'fig2a_tail_ddG_bind_and_tauRAMD')

fig, axs = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
for ax, rows_ in zip(axs, [[('t_a1a_cplx', BLUE, '-', 'TUBA1A tail, Tyr451 (native)'), ('t_a1aYF_cplx', BLUE, '--', 'TUBA1A tail, Y451F'),
                             ('t_a1a_K188A_cplx', GRAY, '-', 'K188A, Tyr451'), ('t_a1aYF_K188A_cplx', GRAY, '--', 'K188A, Y451F')],
                            [('e_a1a_cplx', BLUE, '-', 'TBCE + TUBA1A tail'), ('e_a8_cplx', ORANGE, '--', 'TBCE + TUBA8 tail')]]):
    for sid, c, ls, lab in rows_:
        tt, S, n, med = surv(sid); close(med, ramd_tab[sid]['tau_ns'], 1e-6, f'tauRAMD median {sid}')
        ax.step(np.r_[tt, 12], np.r_[S, S[-1]], where='post', color=c, ls=ls, lw=2.2, label=f'{lab}: τ = {med:.2f} ns')
    ax.axhline(0.5, color=AXIS, lw=1, ls=':'); ax.set_xlim(0, 8); ax.set_xlabel('time under random force (ns)'); ax.legend(fontsize=9.5)
axs[0].set_ylabel('fraction still bound (30 traj.)'); axs[0].set_title('TBCB + TUBA1A tail', loc='left'); axs[1].set_title('TBCE + tails', loc='left')
save(fig, '02_TBCB_tail_clamp', 'fig2b_tauRAMD_other_systems')

# ---------------------------------------------------------------- 3. tail PMF
PM = os.path.join(a.work, 'data', 'results', 'pmf')
fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), sharey=True)
for ax, rec, pairs in ((axs[0], 'TBCE', [('e_a8_cplx_big', ORANGE, 'TUBA8 tail'), ('e_a1a_cplx_big', BLUE, 'TUBA1A tail')]),
                       (axs[1], 'TBCB', [('t_a8_cplx_big', ORANGE, 'TUBA8 tail'), ('t_a1a_cplx_big', BLUE, 'TUBA1A tail')])):
    for sid, c, lab in pairs:
        p = J(os.path.join(PM, sid, 'pmf_results.json'))
        r = np.array(p['r_nm']); W = np.array([np.nan if v is None else v for v in p['W_kcal']], float)
        dW = np.array([np.nan if v is None else v for v in p['dW_kcal']], float)
        Wi = W + KT * np.log(4 * np.pi * r ** 2) - p['bulk_fit_C_kcal']        # interaction PMF, 0 in bulk
        ok_ = p['gate_bulk_reached']; ls = '-' if ok_ else '--'
        ax.plot(r, Wi, ls, color=c, lw=2.4, label=f'{lab}: ΔG° = {fv(p["dG0_kcal"], "+.1f")}' + ('' if ok_ else '  (bulk gate failed)'))
        ax.fill_between(r, Wi - dW, Wi + dW, color=c, alpha=0.12, lw=0)
        if ok_:
            ax.axvspan(r[np.isfinite(W)].max() - 0.6, r[np.isfinite(W)].max(), color=GRID, alpha=0.5, lw=0)
            ax.axvline(p['r_site_boundary_nm'], color=c, lw=1, ls=':')
    ax.axhline(0, color=AXIS, lw=1); ax.set_xlabel('tail – receptor COM distance r (nm)'); ax.set_title(f'{rec} pulling the α tail off', loc='left')
    ax.legend(fontsize=10, loc='upper left' if rec == 'TBCE' else 'lower right')
close(J(os.path.join(PM, 'e_a8_cplx_big', 'pmf_results.json'))['dG0_kcal'], -7.92, 0.006, 'PMF TBCE-TUBA8')
axs[0].set_ylabel('W(r) + kT ln 4πr² − C  (kcal/mol)')
axs[0].text(0.98, 0.04, 'grey band: bulk-fit region\ndotted: site boundary', transform=axs[0].transAxes, fontsize=9, color=INK2, ha='right')
fig.suptitle('Umbrella-sampling PMF (provisional): only TBCE–TUBA8 tail passes the bulk check (K$_d$ ≈ 2.6 µM)', x=0.01, ha='left', fontsize=14, fontweight='bold', y=1.03)
save(fig, '03_tail_PMF_TBCE', 'fig3_tail_PMF')

# ---------------------------------------------------------------- 4. loop 32-45 states in every context
TH = J(I('20261104_2110', 'rest2_state_thresholds.json')); QN, RN, QU = TH['Q12_p05'], TH['rmsd_p95_A'], TH['U_threshold_Q12']
def loopcv(name):
    fs = sorted(glob.glob(I('*', f'loopcv_md310_{name}.json')))
    pref = [f for f in fs if '20261108_final' in f] or fs           # 9M1L: chain B recomputed on 11-08
    return J(pref[-1])
def table_N(name):
    for T in (T8, TM, T7, T6):
        if name in T.get('loop', {}): return T['loop'][name]['N']
CTX = [('free monomer\n200 ns', [('TUBA1A', 'a1a_r0'), ('TUBA1A', 'a1a_r1'), ('TUBA1A', 'a1a_r2'), ('TUBA8', 'a8_r0'), ('TUBA1B', 'a1b_r0'), ('TUBA4A', 'a4a_r0')]),
       ('αβ heterodimer\n(9M1M), 50 ns', [('TUBA1A', 'h_9M1M_Tuba1a_r0'), ('TUBA1A', 'h_9M1M_Tuba1a_r1'), ('TUBA8', 'h_9M1M_Tuba8_r0'), ('TUBA8', 'h_9M1M_Tuba8_r1'),
                                         ('TUBA1B', 'h_9M1M_Tuba1b_r0'), ('TUBA1B', 'h_9M1M_Tuba1b_r1'), ('TUBA4A', 'h_9M1M_Tuba4a_r0'), ('TUBA4A', 'h_9M1M_Tuba4a_r1')]),
       ('TBCE complex\n(9M1L), 50 ns', [('TUBA1A', 'x_9M1L_tbce_a1a_r0'), ('TUBA1A', 'x_9M1L_tbce_a1a_r1'), ('TUBA8', 'x_9M1L_tbce_a8_r0'), ('TUBA8', 'x_9M1L_tbce_a8_r1')]),
       ('TBCE complex\n(AF3 model), 50 ns', [('TUBA1A', 'c_tbce_a1a_m0_r0'), ('TUBA1A', 'c_tbce_a1a_m0_r1'), ('TUBA8', 'c_tbce_a8_m0_r0'), ('TUBA8', 'c_tbce_a8_m0_r1')]),
       ('TBCB complex\n(AF3 model), 50 ns', [('TUBA1A', 'c_tbcb_a1a_m0_r0'), ('TUBA1A', 'c_tbcb_a1a_m0_r1'), ('TUBA1A', 'c_tbcb_a1a_m0_r2'),
                                            ('TUBA8', 'c_tbcb_a8_m0_r0'), ('TUBA8', 'c_tbcb_a8_m0_r1'), ('TUBA8', 'c_tbcb_a8_m0_r2')])]
states = {}
for ctx, items in CTX:
    for iso, name in items:
        d = loopcv(name); q = np.array(d['series']['Q_lam1.2']); rm = np.array(d['series']['rmsd_loop_A']); dt = d['ns_per_frame']
        st = np.where((q >= QN) & (rm <= RN), 0, np.where(q < QU, 2, 1)); after = np.arange(len(q)) * dt >= 10.0
        N = float(np.mean(st[after] == 0)); close(N, table_N(name), 0.003, f'N fraction {name}')
        states[name] = dict(iso=iso, t=np.arange(len(q)) * dt, q=q, st=st, N=N, U=float(np.mean(st[after] == 2)), mirror=d.get('mirror_image', False))
fig, ax = plt.subplots(figsize=(11, 5.2))
yt = []
for k, (ctx, items) in enumerate(CTX):
    y0 = -k
    for j, iso in enumerate(['TUBA1A', 'TUBA8', 'TUBA1B', 'TUBA4A']):
        pts = [states[n] for i_, n in items if i_ == iso]
        if not pts: continue
        off = {'TUBA1A': 0.18, 'TUBA8': -0.18, 'TUBA1B': 0.05, 'TUBA4A': -0.05}[iso]
        for s in pts:
            c = ISO.get(iso, GRAY)
            ax.plot(s['N'], y0 + off, 'o', ms=11 if iso in ISO else 8, color=c, mfc='white' if s['mirror'] else c, mew=2, alpha=1 if iso in ISO else 0.8, zorder=3)
    yt.append((y0, ctx))
ax.set_yticks([q[0] for q in yt]); ax.set_yticklabels([q[1] for q in yt], fontsize=10.5); ax.set_xlim(-0.03, 1.0)
ax.set_xlabel('fraction of time (after 10 ns) with the 32–45 loop in its starting (N) conformation'); ax.grid(axis='y', visible=False)
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, ms=10, label='TUBA1A'), Line2D([], [], marker='o', ls='', color=ORANGE, ms=10, label='TUBA8'),
                   Line2D([], [], marker='o', ls='', color=GRAY, ms=8, label='TUBA1B / TUBA4A'),
                   Line2D([], [], marker='o', ls='', color=ORANGE, mfc='white', mew=2, ms=10, label='mirror-image run (analysed as L)')],
          loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=4, fontsize=10)
ax.set_title('Loop 32–45: TUBA8 never stays native in the heterodimer or the 9M1L TBCE complex, but does in TBCB complexes', loc='left', fontsize=12.5)
save(fig, '04_loop32-45_states', 'fig4a_loop_native_fraction_all_contexts')

# state timelines (N / I / U), every trajectory, one row each
SC = ['#2b2b2a', '#b9b8b0', '#f2f1ec']          # N dark, I mid, U light (ordinal grey ramp; not a series colour)
nrows = sum(len(it) for _, it in CTX) + len(CTX) - 1
fig, ax = plt.subplots(figsize=(13, 0.32 * nrows + 1.6))
row = 0; ylab = []
from matplotlib.colors import ListedColormap
for k, (ctx, items) in enumerate(CTX):
    for iso, name in items:
        s = states[name]
        ax.imshow(s['st'][None, :], extent=[0, s['t'][-1], -row - 0.4, -row + 0.4], aspect='auto', cmap=ListedColormap(SC), vmin=0, vmax=2, interpolation='nearest')
        tag = name.split('_')[-1]
        ylab.append((-row, f'{iso} {tag}' + (' (mirror)' if s['mirror'] else '')))
        ax.text(203, -row, f'N {s["N"]:.2f}', va='center', fontsize=8.5, color=INK2)
        row += 1
    if k < len(CTX) - 1:
        ax.text(-2, -row + 0.05, '', fontsize=1); row += 1
    ax.text(-58, -(row - len(items) - (1 if k < len(CTX) - 1 else 0)) + 0.6, ctx.replace('\n', ' '), fontsize=10, fontweight='bold', color=INK, va='bottom')
ax.set_yticks([q[0] for q in ylab]); ax.set_yticklabels([q[1] for q in ylab], fontsize=8.5)
for t_, lab in zip(ax.get_yticklabels(), ylab):
    t_.set_color(INK)
ax.set_xlim(0, 215); ax.set_ylim(-row + 0.5, 1.6); ax.set_xlabel('time (ns)'); ax.grid(False)
for sp in ('left', 'bottom'): ax.spines[sp].set_color(AXIS)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=SC[0], label=f'N  (Q ≥ {QN:.3f} and loop RMSD ≤ {RN:.2f} Å)'),
                   plt.Rectangle((0, 0), 1, 1, color=SC[1], label='I  (in between)'),
                   plt.Rectangle((0, 0), 1, 1, fc=SC[2], ec=AXIS, label=f'U  (Q < {QU:.3f})')], loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=10)
save(fig, '04_loop32-45_states', 'fig4b_loop_state_timelines')

# Q(t) traces: TUBA1A vs TUBA8 per context, 1-ns running mean
fig, axs = plt.subplots(len(CTX), 2, figsize=(13, 2.1 * len(CTX)), sharey=True)
for k, (ctx, items) in enumerate(CTX):
    for j, iso in enumerate(['TUBA1A', 'TUBA8']):
        ax = axs[k, j]
        for i_, n in items:
            if i_ != iso: continue
            s = states[n]; w = max(1, int(round(1.0 / (s['t'][1] - s['t'][0]))))
            qm = rmean(s['q'], w)
            ax.plot(s['t'], qm, color=ISO[iso], lw=1.4, alpha=0.9)
            ax.text(s['t'][-1] * 1.01, qm[-w], n.split('_')[-1], fontsize=8, color=INK2, va='center')
        ax.axhline(QN, color=INK2, lw=0.8, ls='--'); ax.axhline(QU, color=INK2, lw=0.8, ls=':')
        ax.set_ylim(0, 1.0); ax.set_xlim(0, (200 if k == 0 else 50) * 1.08)
        if k == 0: ax.set_title(iso, loc='left', color=INK)
        if j == 0: ax.set_ylabel({0: 'free monomer', 1: 'αβ heterodimer', 2: 'TBCE (9M1L)', 3: 'TBCE (AF3)', 4: 'TBCB (AF3)'}[k], fontsize=10)
axs[-1, 0].set_xlabel('time (ns)'); axs[-1, 1].set_xlabel('time (ns)')
fig.suptitle(f'Loop 32–45 native-contact fraction Q(t), 1-ns running mean (dashed: Q part of the N test, {QN:.2f}; dotted: U threshold {QU:.2f})\n'
             f'N also needs loop RMSD ≤ {RN:.2f} Å, so Q above the dashed line alone is not N (state calls: fig4b)',
             x=0.01, ha='left', fontsize=13, fontweight='bold', y=1.01)
fig.tight_layout()
save(fig, '04_loop32-45_states', 'fig4c_loop_Q_timeseries')

# ---------------------------------------------------------------- 5. loop residue energetics
RES = T6['residue']; RM = J(I('20261110_m', 'residue_a1a_10-200ns.json'))['per_residue']
def sysv(sysname, pos, q):
    p = RES[f'{sysname}.json']['positions'][str(pos)]
    return {'E': (p['E_local'], p['E_se']), 'LJ': (p['E_lj'], p['E_lj_se']), 'SASA': (p['rel_sasa'], p['SASA_se'])}[q], p['resname']
def a1a_r12(pos, q):
    p = RM[str(pos)]
    if q == 'E': return p['E_coul_kcal'][0] + p['E_lj_kcal'][0], np.hypot(p['E_coul_kcal'][1], p['E_lj_kcal'][1])
    if q == 'LJ': return tuple(p['E_lj_kcal'])
    return tuple(p['rel_sasa_mean'])
PANELS = [(39, 'E', 'residue 39\nlocal energy (kcal/mol)'), (42, 'E', 'residue 42\nlocal energy (kcal/mol)'), (45, 'E', 'residue 45\nlocal energy (kcal/mol)'),
          (43, 'SASA', 'residue 43\nrelative SASA (%)'), (114, 'LJ', 'residue 114\nLJ energy (kcal/mol)')]
close(sysv('a8', 39, 'E')[0][0], -29.6, 0.06, 'E39 TUBA8'); close(sysv('a1a', 39, 'E')[0][0], 118.6, 0.06, 'E39 TUBA1A r0')
close(a1a_r12(39, 'E')[0], 50.0, 0.06, 'E39 TUBA1A r1/r2'); close(sysv('a8', 114, 'LJ')[0][0], -10.4, 0.06, 'LJ114 TUBA8')
fig, axs = plt.subplots(1, len(PANELS), figsize=(16, 4.4))
for ax, (pos, q, lab) in zip(axs, PANELS):
    pts = [('TUBA8', 0, *sysv('a8', pos, q), ORANGE), ('TUBA1A r0', 1, *sysv('a1a', pos, q), BLUE)]
    v12 = a1a_r12(pos, q)
    for name, x, (v, e), rn, c in pts:
        ax.errorbar(x, v, yerr=e, fmt='o', ms=10, color=c, elinewidth=2); ax.text(x + 0.12, v, f'{rn}\n{fv(v, ".0f")}' if q != 'SASA' else f'{rn}\n{v:.0f}%', fontsize=9, va='center', color=INK)
    ax.errorbar(1.45, v12[0], yerr=v12[1], fmt='s', ms=9, color=BLUE, mfc='white', mew=2, elinewidth=2)
    ax.text(1.57, v12[0], fv(v12[0], '.0f'), fontsize=9, va='center', color=INK)
    for x, sname in ((2.3, 'a1b'), (2.75, 'a4a')):
        (v, e), rn = sysv(sname, pos, q); ax.errorbar(x, v, yerr=e, fmt='o', ms=8, color=GRAY, elinewidth=2)
    ax.set_xticks([0, 1, 1.45, 2.3, 2.75]); ax.set_xticklabels(['A8', '1A\nr0', '1A\nr1,2', '1B', '4A'], fontsize=9.5)
    ax.set_xlim(-0.5, 3.2); ax.set_title(lab, loc='left', fontsize=10.5); ax.grid(axis='x', visible=False)
    if q != 'SASA': ax.axhline(0, color=AXIS, lw=1)
fig.suptitle('TUBA8 loop substitutions rewire the local energetics (10–200 ns; A8 = TUBA8, 1A/1B/4A = TUBA1A/1B/4A; bars: block SE, or SD over r1/r2)',
             x=0.01, ha='left', fontsize=13, fontweight='bold', y=1.04)
save(fig, '05_loop_local_energetics', 'fig5_loop_residue_energetics')

# ---------------------------------------------------------------- 7. MM/GBSA
mm = T8['mmgbsa']
G = [('TBCB (AF3 model)', [('TUBA1A', 'c_tbcb_a1a_m0'), ('TUBA8', 'c_tbcb_a8_m0')]),
     ('TBCE (AF3 model)', [('TUBA1A', 'c_tbce_a1a_m0'), ('TUBA8', 'c_tbce_a8_m0')]),
     ('TBCE (9M1L)', [('TUBA1A', 'x_9M1L_tbce_a1a'), ('TUBA8', 'x_9M1L_tbce_a8')]),
     ('αβ heterodimer (9M1M)', [('TUBA1A', 'h_9M1M_Tuba1a'), ('TUBA8', 'h_9M1M_Tuba8'), ('TUBA1B', 'h_9M1M_Tuba1b'), ('TUBA4A', 'h_9M1M_Tuba4a')])]
close(mm['c_tbcb_a1a_m0']['GB']['mean'], -153, 0.5, 'MMGBSA TBCB-TUBA1A'); close(mm['c_tbcb_a8_m0']['GB']['mean'], -176, 0.5, 'MMGBSA TBCB-TUBA8')
fig, ax = plt.subplots(figsize=(10, 4.6)); xt = []; x = 0
for grp, items in G:
    x0 = x
    for iso, sid in items:
        gb = mm[sid]['GB']; c = ISO.get(iso, GRAY)
        ax.errorbar(x, gb['mean'], yerr=gb['sd_replicas'], fmt='o', ms=10, color=c, elinewidth=2)
        ax.text(x, gb['mean'] - gb['sd_replicas'] - 6, f'{fv(gb["mean"], ".0f")}\n(n={gb["n_replicas"]})', ha='center', va='top', fontsize=8.5, color=INK2)
        x += 1
    xt.append(((x0 + x - 1) / 2, grp)); x += 0.8
ax.set_xticks([t[0] for t in xt]); ax.set_xticklabels([t[1] for t in xt], fontsize=10); ax.grid(axis='x', visible=False)
ax.set_ylabel('MM/GBSA score (kcal/mol), mean ± SD over replicas'); ax.set_ylim(-245, 0)
ax.legend(handles=[Line2D([], [], marker='o', ls='', color=BLUE, ms=9, label='TUBA1A'), Line2D([], [], marker='o', ls='', color=ORANGE, ms=9, label='TUBA8'),
                   Line2D([], [], marker='o', ls='', color=GRAY, ms=9, label='TUBA1B / TUBA4A')], loc='lower right', fontsize=10, ncol=3)
ax.set_title('MM/GBSA (ranking score, not an affinity): no isoform difference beyond replica spread', loc='left', fontsize=12.5)
save(fig, '07_MMGBSA_ranking', 'fig7_mmgbsa')

# ---------------------------------------------------------------- figures from the a31 clip series
SER = {}
for d in (os.path.join(a.work, 'out_sherlock'), os.path.join(a.work, 'out_marlowe')):
    for f in glob.glob(os.path.join(d, '*_series.json')):
        SER[os.path.basename(f)[:-12]] = J(f)
print('series available:', sorted(SER))
def rep1(label):
    r = SER[label]['replicas']; return list(r.values())[0]

if all(k in SER for k in ('free_TUBA8_r0', 'free_TUBA1A_r0', 'free_TUBA1A_I114S_r0')):
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.4), gridspec_kw={'width_ratios': [1.6, 1]})
    ax = axs[0]
    for lab, c, name in (('free_TUBA8_r0', ORANGE, 'TUBA8 (Ser114)'), ('free_TUBA1A_I114S_r0', BLUE, 'TUBA1A-I114S (Ser114)')):
        r = rep1(lab); d = np.array(r['res114']['d_OG114_O110_A']); t = np.arange(len(d)) * r['ns_per_frame']
        w = int(round(1 / r['ns_per_frame'])); dm = rmean(d, w)
        ax.plot(t, d, color=c, lw=0.3, alpha=0.25); ax.plot(t, dm, color=c, lw=1.8, label=f'{name}: {np.mean(d[t >= 10] <= 3.5) * 100:.0f}% of frames ≤ 3.5 Å')
    ax.axhline(3.5, color=INK2, lw=1, ls='--'); ax.set_ylim(2, 9); ax.set_xlabel('time (ns)'); ax.set_ylabel('d(Ser114 OG – Val110 O)  (Å)')
    ax.set_title('Ser114 hydroxyl ↔ Val110 backbone O (dashed: 3.5 Å H-bond cut-off)', loc='left'); ax.legend(fontsize=10, loc='upper right')
    ax = axs[1]; groups = [('TUBA8\nSer114', ['free_TUBA8_r0'], ORANGE), ('TUBA1A-I114S\nSer114', ['free_TUBA1A_I114S_r0'], BLUE),
                           ('TUBA1A\nIle114', ['free_TUBA1A_r0', 'free_TUBA1A_r1', 'free_TUBA1A_r2'], BLUE), ('TUBA1B / 4A\nIle114', ['free_TUBA1B_r0', 'free_TUBA4A_r0'], GRAY)]
    for x, (lab, labs, c) in enumerate(groups):
        for k, l in enumerate([l for l in labs if l in SER]):
            r = rep1(l); n = np.array(r['res114']['n_contacts_sc_4.5A']); t = np.arange(len(n)) * r['ns_per_frame']; n = n[t >= 10]
            xx = x + (k - (len(labs) - 1) / 2) * 0.18
            ax.plot([xx, xx], np.percentile(n, [25, 75]), color=c, lw=5, solid_capstyle='round', alpha=0.6)
            ax.plot(xx, np.median(n), 'o', color=c, ms=9, mfc='white' if 'Ile' in lab else c, mew=2)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g_[0] for g_ in groups], fontsize=9.5); ax.grid(axis='x', visible=False)
    ax.set_ylabel('heavy atoms within 4.5 Å of the 114 side chain'); ax.set_title('Packing around residue 114 (median, IQR)', loc='left')
    save(fig, '01_Ser114_folding_cost', 'fig1d_Ser114_hbond_and_packing_MD')

tb = [l for l in ('TBCB_TUBA8_r1', 'TBCB_TUBA8_r2', 'TBCB_TUBA1A_r0', 'TBCB_TUBA1A_r1', 'TBCB_TUBA1A_r2') if l in SER]
if len(tb) == 5:
    fig, axs = plt.subplots(2, 1, figsize=(13, 6.4), sharex=True)
    for ax, iso, labs in ((axs[0], 'TUBA8', tb[:2]), (axs[1], 'TUBA1A', tb[2:])):
        for k, l in enumerate(labs):
            r = rep1(l); gg = r['groove']; t = np.arange(len(gg['d_ring_PHE187_A'])) * r['ns_per_frame']
            dmin = np.minimum(np.array(gg['d_ring_PHE187_A']), np.array(gg['d_ring_TRP192_A']))
            w = int(round(0.5 / r['ns_per_frame'])); dm = rmean(dmin, w)
            ax.plot(t, dm, color=ISO[iso], lw=1.8, alpha=[1, 0.7, 0.45][k], label=f'{l.split("_")[-1]}: {np.mean(dmin <= 6.5) * 100:.0f}% ≤ 6.5 Å')
            if 'd_OH451_NZ188_A' in gg:
                oh = rmean(gg['d_OH451_NZ188_A'], w)
                ax.plot(t, oh, color=GRAY, lw=1, ls='--', alpha=[1, 0.7, 0.45][k])
        ax.axhline(6.5, color=INK2, lw=1, ls=':'); ax.set_ylim(3, 25); ax.set_ylabel('distance (Å)')
        ax.set_title(f'TBCB + {iso}: tail {gg["tail_residue"]} ring to the nearer of F187 / W192 rings' + ('   (grey dashed: Tyr451 OH – K188 NZ)' if iso == 'TUBA1A' else ''), loc='left', fontsize=11.5)
        ax.legend(fontsize=9.5, loc='upper right', ncol=3)
    axs[1].set_xlabel('time (ns), 0.5-ns running mean')
    fig.suptitle('The CAP-Gly groove (F187 / W192) holds the tail-end aromatic ring', x=0.01, ha='left', fontsize=14, fontweight='bold', y=1.01)
    fig.tight_layout(); save(fig, '02_TBCB_tail_clamp', 'fig2c_groove_ring_distances_MD')

U8 = [l for l in ('md450_TUBA1A', 'md450_TUBA8', 'md450_TUBA1A_I114S', 'md450_TUBA1B', 'md450_TUBA4A', 'md450_TUBA8_S114I') if l in SER]
if len(U8) == 6:
    GRP = ['loop32_45', 'unit1', 'unit2', 'unit3', 'unit0', 'seg110_120', 'unit5', 'unit4']
    GL = {'loop32_45': 'loop 32–45', 'seg110_120': 'seg 110–120', **{f'unit{i}': f'unit {i}' for i in range(6)}}
    ref18 = {'md450_TUBA1A': T6['unfold']['a1a'], 'md450_TUBA8': T6['unfold']['a8'], 'md450_TUBA1A_I114S': TM['unfold']['a1a_I114S'],
             'md450_TUBA1B': TM['unfold']['a1b'], 'md450_TUBA4A': TM['unfold']['a4a'], 'md450_TUBA8_S114I': TM['unfold']['a8_S114I']}
    fig, axs = plt.subplots(2, 3, figsize=(15, 6.4), sharex=True, sharey=True)
    thalf = {}
    for ax, l in zip(axs.flat, U8):
        reps = SER[l]['replicas']; Qm = []
        for rep, r in reps.items():
            dt = r['ns_per_frame']; w = int(round(1.0 / dt)); rowq = []
            for gname in GRP:
                q = np.array(r['units'][gname]); rm = np.convolve(q, np.ones(w) / w, mode='valid')
                below = np.where(rm < 0.5)[0]; th = float((below[0] + w / 2) * dt) if len(below) else None
                th18 = ref18[l]['replicas'][rep][gname]['t_half_ns']
                assert (th is None and th18 is None) or abs(th - th18) < 0.02, (l, rep, gname, th, th18)
                thalf[(l, rep, gname)] = th
                rowq.append(rmean(q, w))
            Qm.append(rowq)
        Qm = np.mean(Qm, axis=0)
        im = ax.imshow(Qm, aspect='auto', extent=[0, 20, len(GRP) - 0.5, -0.5], cmap='Blues', vmin=0.2, vmax=0.9, interpolation='nearest')
        ax.set_title(l.replace('md450_', '').replace('_', '-') + ('  (Ser114)' if l in ('md450_TUBA8', 'md450_TUBA1A_I114S') else '  (Ile114)'), loc='left', fontsize=11)
        ax.set_yticks(range(len(GRP))); ax.set_yticklabels([GL[g_] for g_ in GRP], fontsize=9); ax.grid(False)
    for ax in axs[-1]: ax.set_xlabel('time at 450 K (ns)')
    cb = fig.colorbar(im, ax=axs, shrink=0.8, pad=0.01); cb.set_label('native-contact fraction Q (mean of 3 replicas, 1-ns mean)')
    fig.suptitle('450 K unfolding: the 32–45 loop goes first, units 4 and 5 last, in all six systems', x=0.01, ha='left', fontsize=14, fontweight='bold', y=1.0)
    save(fig, '06_450K_unfolding_order', 'fig6a_450K_unit_Q_heatmaps')
    fig, ax = plt.subplots(figsize=(11, 4.8))
    for j, gname in enumerate(GRP):
        for l in U8:
            ser = l in ('md450_TUBA8', 'md450_TUBA1A_I114S')
            for k, rep in enumerate(sorted(SER[l]['replicas'])):
                th = thalf[(l, rep, gname)]; jit = (U8.index(l) * 3 + k - 8.5) * 0.035
                if th is None: ax.plot(21.2, j + jit, '>', color=ORANGE if ser else BLUE, mfc='white', ms=6, mew=1.5, alpha=0.8)
                else: ax.plot(th, j + jit, 'o', color=ORANGE if ser else BLUE, ms=6, alpha=0.85)
    ax.axvline(20, color=AXIS, lw=1); ax.text(20.3, -0.8, 'not reached\nin 20 ns', fontsize=9, color=INK2)
    ax.set_yticks(range(len(GRP))); ax.set_yticklabels([GL[g_] for g_ in GRP]); ax.invert_yaxis(); ax.set_xlim(0, 22.5); ax.grid(axis='y', visible=False)
    ax.set_xlabel('t½: first time the 1-ns mean of Q falls below 0.5 (ns)')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=ORANGE, ms=8, label='Ser at 114 (TUBA8, TUBA1A-I114S)'),
                       Line2D([], [], marker='o', ls='', color=BLUE, ms=8, label='Ile at 114 (TUBA1A, 1B, 4A, TUBA8-S114I)'),
                       Line2D([], [], marker='>', ls='', color=GRAY, mfc='white', ms=8, label='censored at 20 ns')], loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=9.5)
    ax.set_title('Per-replica t½ (18 replicas): same order everywhere; 110–120 not specific to Ser114', loc='left', fontsize=12.5)
    save(fig, '06_450K_unfolding_order', 'fig6b_450K_t_half_per_replica')

print('\n'.join(made))
