#!/usr/bin/env python3
"""Number tables for the interim results report, read only from complete analysis JSON files (never from logs).
One generator for every number in the report: FEP legs and ddG (a09 + a27), loop 32-45 states (a26 series with the
calibrated thresholds), complex interfaces (a25), residue level (a13), 450 K unfolding (a18), REST2 progress (a10),
docking controls (HADDOCK caprieval JSON), tail binding ddG and tauRAMD (a22), US-phi43 (a17), MSM (a11).
Usage: compile_interim_tables.py --interim <dir> [--iface_dir <dir>] [--dock_dir <dir>] --thresholds <json>
                                 [--fep_extra LABEL=<fep_interim.json>] --out <prefix>   (writes <prefix>.md and .json)
"""
import argparse, glob, json, math, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--interim', required=True)
ap.add_argument('--iface_dir', default=None, help='directory with iface_*.json (default: --interim)')
ap.add_argument('--dock_dir', default=None)
ap.add_argument('--thresholds', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--conv_extra', action='append', default=[], help='LABEL=fep_convergence.json of an earlier pass (repeatable)')
ap.add_argument('--residue_extra', action='append', default=[], help='LABEL=residue_*.json of an earlier pass (repeatable)')
ap.add_argument('--fep_extra', action='append', default=[], help='LABEL=fep_interim.json of an earlier pass (repeatable)')
ap.add_argument('--unfold_extra', action='append', default=[], help='SID=a18 json (e.g. results/unfold/a8.json) (repeatable)')
ap.add_argument('--rest2_extra', action='append', default=[], help='SID=rest2_results.json (repeatable)')
ap.add_argument('--msm', action='append', default=[], help='LABEL=msm_results.json (repeatable)')
ap.add_argument('--pmf_dir', default=None, help='results dir with pmf/*/ and pmfdom/*/pmf_results.json')
ap.add_argument('--mmgbsa_dir', default=None, help='directory with a12b <system>.json')
ap.add_argument('--binding_dir', default=None, help='directory with tail_rbfe.json / ramd_results.json (default: <interim>/binding)')
a = ap.parse_args()
I = a.interim
OUT = {}
md = []


def f(x, n=2):
    return '—' if x is None or (isinstance(x, float) and math.isnan(x)) else f'{x:.{n}f}'


def load(p):
    with open(p) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- FEP
def ext(R):
    """A replicate whose checkpoint has more iterations than its status target is being extended: progress value."""
    st = R.get('status') or {}
    tgt = st.get('target')
    return bool(tgt) and R['n_iter'] > tgt


def fep_section(fep, label, key):
    md.append(f'## FEP：{label}（HREX + MBAR；kcal/mol，310 K）\n')
    md.append('| 臂 | 副本 | 迭代 | 状态 | ΔG ± SE | 前半 / 后半 | 半程差 | 最小相邻重叠 | 往返次数 |')
    md.append('|---|---|---|---|---|---|---|---|---|')
    O = OUT.setdefault('fep', {}).setdefault(key, {})
    for leg, L in fep['legs'].items():
        O[leg] = {'summary': L['summary'], 'replicates': {}}
        for r, R in sorted(L['replicates'].items()):
            st = R.get('status') or {}
            done = bool(st.get('done'))
            tgt = st.get('target')
            state = '延长中（进度值）' if ext(R) else ('完成' if done else '进行中')
            h1, h2 = R.get('dG_first_half'), R.get('dG_second_half')
            hs = abs(h1 - h2) if (h1 is not None and h2 is not None) else None
            rt = R.get('round_trips')
            rt = int(np.sum(rt)) if isinstance(rt, list) else rt
            O[leg]['replicates'][r] = dict(n_iter=R['n_iter'], target=tgt, state=state, dG=R['dG_kcal'],
                                           se=R['se_kcal'], half1=h1, half2=h2, half_split=hs,
                                           min_overlap=R['min_adjacent_overlap'], round_trips=rt)
            md.append(f"| {leg} | {r} | {R['n_iter']}{'（原目标 ' + str(tgt) + '）' if tgt else ''} | {state} | "
                      f"{f(R['dG_kcal'], 3)} ± {f(R['se_kcal'], 3)} | {f(h1, 3)} / {f(h2, 3)} | {f(hs, 3)} | "
                      f"{f(R['min_adjacent_overlap'], 3)} | {rt} |")
    md.append('')
    md.append('| 臂 | 计入副本数 | 均值 | SD | SEM | 最大半程差 | 半程门 | 含延长中副本 |')
    md.append('|---|---|---|---|---|---|---|---|')
    for leg, L in fep['legs'].items():
        s = L['summary']
        g = s.get('gate_half_split')
        e = any(ext(R) for R in L['replicates'].values())
        md.append(f"| {leg} | {s['n_passing']}/{s['n_replicates']} | {f(s['mean_dG_kcal'], 3)} | {f(s['sd_dG_kcal'], 3)} | "
                  f"{f(s['sem_dG_kcal'], 3)} | {f(s.get('max_half_split_kcal'), 3)} | "
                  f"{'—' if g is None else ('通过' if g else '未通过')} | {'是' if e else '否'} |")
    md.append('')

    def leg_mean(leg):
        L = fep['legs'].get(leg)
        if not L:
            return None, None, 0
        v = [R['dG_kcal'] for R in L['replicates'].values()]
        if not v:
            return None, None, 0
        return float(np.mean(v)), (float(np.std(v, ddof=1)) if len(v) > 1 else None), len(v)

    D_ = OUT.setdefault('ddG', {}).setdefault(key, {})
    md.append('| 突变 | ΔΔG_fold = ΔG_折叠 − ΔG_非折叠 | 不确定度 | 依据 |')
    md.append('|---|---|---|---|')
    for k, D in fep.get('ddG_fold', {}).items():
        legs = [x.strip() for x in k.replace('[ref', '|').replace(']', '').split('|')]
        fl = legs[0] + '_folded'
        e = any(ext(R) for lg in (fl, legs[-1]) for R in (fep['legs'].get(lg) or {}).get('replicates', {}).values())
        nf, nu = D['folded']['n_passing'], D['unfolded']['n_passing']
        D_[k] = dict(D, includes_extending=e)
        unc = (f"±{f(D['sd_kcal'], 2)}（两臂副本 SD 合并）" if min(nf, nu) >= 2 else '—（有一臂只有 1 个完成副本）')
        md.append(f"| {k} | {f(D['ddG_fold_kcal'], 2)}{'（进度值）' if e else ''} | {unc} | "
                  f"折叠 {nf} + 非折叠 {nu} 个计入副本 |")
    # provisional pair with every replicate that has data (incomplete ones included, clearly labelled)
    for mut, fl, ul in [('S114I（TUBA8，含未完成副本，暂定）', 'a8_114ILE_folded', 'a8_114ILE_unfolded')]:
        mf, sf, nf = leg_mean(fl)
        mu, su, nu = leg_mean(ul)
        if mf is not None and mu is not None:
            d = mf - mu
            unc = math.sqrt((sf or 0) ** 2 + (su or 0) ** 2) if (sf or su) else None
            D_[mut] = dict(ddG_fold_kcal=d, sd_kcal=unc, n_folded=nf, n_unfolded=nu, provisional=True)
            md.append(f"| {mut} | {f(d, 2)} | {('±' + f(unc, 2)) if unc else '—'}（副本 SD 合并） | "
                      f"折叠 {nf} 个、非折叠 {nu} 个副本（含未完成） |")
    md.append('')


fep = load(os.path.join(I, 'fep_interim.json'))
conv_p = os.path.join(I, 'fep_convergence.json')
conv = load(conv_p) if os.path.exists(conv_p) else None
fep_section(fep, f'本轮（{os.path.basename(os.path.normpath(I))}）', 'this_pass')
for spec in a.fep_extra:
    lab, path = spec.split('=', 1)
    fep_section(load(path), lab, lab)

# a27: time convergence and the path decomposition of the ddG
def a27_section(conv, label, key):
    C = conv['legs']
    md.append(f'### a27 收敛诊断：{label}（只读检查点）\n')
    md.append('| 臂 | 副本 | 4 个连续块 ΔG | 末 20/40/60/80/100 % 生产段 ΔG | pymbar 平衡点 t_eq / N_eff |')
    md.append('|---|---|---|---|---|')
    A27 = OUT.setdefault('a27', {}).setdefault(key, {})
    for leg in sorted(C):
        for r, R in sorted(C[leg].items()):
            if not isinstance(R, dict) or 'blocks' not in R:
                continue
            bl = ' / '.join(f(b['dG'], 2) for b in R['blocks'])
            rv = ' / '.join(f(R['reverse'][k]['dG'], 2) for k in ['0.2', '0.4', '0.6', '0.8', '1.0'] if k in (R.get('reverse') or {}))
            pe = R.get('pymbar_equilibration') or {}
            md.append(f"| {leg} | {r} | {bl} | {rv} | {pe.get('t_eq')} / {f(pe.get('N_eff'), 0)} |")
            A27.setdefault(leg, {})[r] = dict(blocks=[b['dG'] for b in R['blocks']],
                                                      reverse={k: v['dG'] for k, v in (R.get('reverse') or {}).items()},
                                                      steps={k: v['full'] for k, v in (R.get('steps') or {}).items()})
    md.append('')
    steps = []
    for leg in sorted(C):
        for r, R in sorted(C[leg].items()):
            for k in (R.get('steps') or {}):
                if k not in steps:
                    steps.append(k)
    if steps:
        md.append('每个副本沿 λ 路径的分步 ΔG（全生产段）：\n')
        md.append('| 臂 | 副本 | ' + ' | '.join(steps) + ' |')
        md.append('|---|---|' + '---|' * len(steps))
        for leg in sorted(C):
            for r, R in sorted(C[leg].items()):
                st = R.get('steps') or {}
                if st:
                    md.append(f'| {leg} | {r} | ' + ' | '.join(f(st[k]['full'], 2) if k in st else '—' for k in steps) + ' |')
        md.append('')
    for mut, fl, ul in [('I114S', 'a1a_114SER_folded', 'a1a_114SER_unfolded'),
                        ('S114I（暂定）', 'a8_114ILE_folded', 'a8_114ILE_unfolded')]:
        if fl in C and ul in C:
            def stepmean(leg):
                acc = {}
                for r, R in C[leg].items():
                    for k, v in (R.get('steps') or {}).items():
                        acc.setdefault(k, []).append(v['full'])
                return {k: float(np.mean(v)) for k, v in acc.items()}
            sf_, su_ = stepmean(fl), stepmean(ul)
            if sf_ and su_:
                md.append(f'{mut} 沿 λ 路径分解（折叠减非折叠，副本平均；依赖路径的拆分，不是唯一的热力学分解）：')
                rows = []
                for k in sf_:
                    if k in su_:
                        rows.append((k, sf_[k] - su_[k]))
                md.append('；'.join(f'{k} {f(v, 2)}' for k, v in rows) + f'；合计 {f(sum(v for _, v in rows), 2)}\n')
                OUT.setdefault('path_decomposition', {}).setdefault(key, {})[mut] = dict(rows)
                # ddG from the last x % of each leg
                rev = {}
                for x in ['0.2', '0.4', '0.6', '0.8', '1.0']:
                    vf = [R['reverse'][x]['dG'] for R in C[fl].values() if x in (R.get('reverse') or {})]
                    vu = [R['reverse'][x]['dG'] for R in C[ul].values() if x in (R.get('reverse') or {})]
                    if vf and vu:
                        rev[x] = float(np.mean(vf) - np.mean(vu))
                md.append(f'{mut} 用每条臂末 20/40/60/80/100 % 生产段：' +
                          ' / '.join(f(rev[x], 2) for x in rev) + '\n')
                OUT.setdefault('ddG_reverse_cumulative', {}).setdefault(key, {})[mut] = rev

if conv:
    a27_section(conv, f'本轮（{os.path.basename(os.path.normpath(I))}）', 'this_pass')
for spec in a.conv_extra:
    lab, path = spec.split('=', 1)
    a27_section(load(path), lab, lab)

# ---------------------------------------------------------------- loop 32-45 states
thr = load(a.thresholds)
QN, RN, QU = thr['Q12_p05'], thr['rmsd_p95_A'], thr['U_threshold_Q12']
md.append('## 环 32–45 状态（310 K 无偏 MD；N: Q(λ=1.2) ≥ %.3f 且 RMSD ≤ %.2f Å；U: Q(λ=1.2) < %.3f；其余为 I）\n'
          % (QN, RN, QU))
md.append('| 体系 | 轨迹 | 长度 ns | 10 ns 后 N / I / U | 每 20 ns 窗 N 比例 | 每 20 ns 窗 U 比例 | 末 10 ns 环 RMSD Å | 末 10 ns Q(λ=1.8) | φ43>0 比例 |')
md.append('|---|---|---|---|---|---|---|---|---|')
OUT['loop'] = {}
for p in sorted(glob.glob(os.path.join(I, 'loopcv_md310_*.json'))):
    d = load(p)
    sid = os.path.basename(p)[len('loopcv_md310_'):-5]
    s = d['series']
    q = np.array(s['Q_lam1.2'], float); rm = np.array(s['rmsd_loop_A'], float)
    q8 = np.array(s['Q_lam1.8'], float); phi = np.array(s['phi43'], float)
    dt = d['ns_per_frame']; t = (np.arange(len(q)) + 1) * dt
    N = (q >= QN) & (rm <= RN); U = q < QU
    m = t > 10.0
    win = []
    for w0 in np.arange(0, t[-1], 20.0):
        k = (t > w0) & (t <= w0 + 20.0)
        if k.sum() > 0:
            win.append((float(w0), float(N[k].mean()), float(U[k].mean()), int(k.sum())))
    last = t > t[-1] - 10.0
    rec = dict(traj='/'.join(d['dir'].split('/')[-2:]), ns=float(t[-1]), n_frames=int(len(q)),
               mirror_image=d.get('mirror_image'), res43=d.get('res43'),
               N=float(N[m].mean()), U=float(U[m].mean()), I=float(1 - N[m].mean() - U[m].mean()),
               windows=win, last10_rmsd=float(rm[last].mean()), last10_Q18=float(q8[last].mean()),
               p_phi43_pos=float((phi > 0).mean()))
    OUT['loop'][sid] = rec
    md.append(f"| {sid} | {rec['traj']}{'（镜像轨迹，φ/ψ 已换回 L 系）' if rec['mirror_image'] else ''} | {f(rec['ns'], 1)} | "
              f"{f(rec['N'])} / {f(rec['I'])} / {f(rec['U'])} | {' '.join(f(w[1]) for w in win)} | "
              f"{' '.join(f(w[2]) for w in win)} | {f(rec['last10_rmsd'])} | {f(rec['last10_Q18'])} | "
              f"{f(rec['p_phi43_pos'])}（{rec['res43']}43） |")
md.append('')

# ---------------------------------------------------------------- interfaces
ID = a.iface_dir or I
OUT['iface'] = {}
ifs = sorted(glob.glob(os.path.join(ID, 'iface_*.json')))
if ifs:
    md.append('## 复合物界面（a25，310 K，四个等长块）\n')
    md.append('| 体系 | 轨迹 | 重原子接触数（块） | Q_if（块） | 界面 CA RMSD（叠合 α 后，块）Å | 伙伴骨架 RMSD（块）Å | 最小镜像差 |')
    md.append('|---|---|---|---|---|---|---|')
    for p in ifs:
        d = load(p)
        sid = d['system']
        for r, R in d['replicas'].items():
            OUT['iface'][f'{sid}/{r}'] = {k: R[k] for k in R if k in ('heavy_contacts_blocks', 'Q_if_blocks',
                                          'iface_ca_rmsd_on_A_blocks_A', 'rmsd_bb_B_blocks_A', 'trajectory_ns',
                                          'interchain_salt_bridges_second_half', 'focus', 'minimum_image_max_difference')}
            md.append(f"| {sid} | {r}，{f(R['trajectory_ns'], 1)} ns | {' / '.join(f(x, 0) for x in R['heavy_contacts_blocks'])} | "
                      f"{' / '.join(f(x) for x in R['Q_if_blocks'])} | {' / '.join(f(x, 1) for x in R['iface_ca_rmsd_on_A_blocks_A'])} | "
                      f"{' / '.join(f(x, 1) for x in R['rmsd_bb_B_blocks_A'])} | {R.get('minimum_image_max_difference')} |")
    md.append('')
    for p in ifs:
        d = load(p)
        for r, R in d['replicas'].items():
            fo = R.get('focus') or {}
            if fo:
                md.append(f"{d['system']} {r} 关注残基（前半 → 后半占有率）：")
                for k, v in fo.items():
                    fh = {x[0]: x[1] for x in v['first_half']}; sh = {x[0]: x[1] for x in v['second_half']}
                    partners = sorted(set(fh) | set(sh))
                    md.append(f"- {k} {v['resname']}: " + '，'.join(f"{pp} {f(fh.get(pp), 2)}→{f(sh.get(pp), 2)}" for pp in partners))
                md.append('')
            sb = R.get('interchain_salt_bridges_second_half') or []
            if sb:
                md.append(f"{d['system']} {r} 后半段盐桥（占有率 ≥ 0.3）：" +
                          '，'.join(f'{x[0]}–{x[1]} {f(x[2], 2)}' for x in sb if x[2] >= 0.3) + '\n')

# ---------------------------------------------------------------- residue level
res = sorted(glob.glob(os.path.join(I, 'residue_*ns.json'))) + [x.split('=', 1)[1] for x in a.residue_extra]
POS = ['32', '39', '40', '41', '42', '43', '45', '46', '110', '114']
if res:
    md.append('## 残基层（a13；局部相互作用能 = 该残基与 1.2 nm 内蛋白 + GTP·Mg 的库仑（介电 1）+ LJ，是势能项不是自由能；'
              'SASA 为相对 SASA 百分数）\n')
    md.append('| 文件 | 窗口 ns | 周期镜像检查（nm） | ' + ' | '.join(f'{p}' for p in POS) + ' |')
    md.append('|---|---|---|' + '---|' * len(POS))
    OUT['residue'] = {}
    for p in res:
        d = load(p)
        pr = d['per_residue']
        w = d.get('window_ns_frames_per_rep', [[None, None, None]])[0]
        img = d.get('image_check_max_nm_per_rep')
        cells = []
        rec = {}
        for k in POS:
            if k in pr:
                x = pr[k]
                se = x['block_se_per_rep'][0]
                ec, el = x['E_coul_kcal'][0], x['E_lj_kcal'][0]
                e = (ec or 0) + (el or 0)
                ese = math.hypot(se.get('E_coul_block_se_kcal', 0) or 0, se.get('E_lj_block_se_kcal', 0) or 0)
                rec[k] = dict(resname=x['resname'], E_local=e, E_se=ese, E_coul=ec, E_lj=el,
                              E_lj_se=se.get('E_lj_block_se_kcal'), rel_sasa=x['rel_sasa_mean'][0], SASA_se=se.get('rel_sasa_block_se'),
                              hydrophobic_partners=(x.get('hydrophobic_partners_mean') or [None])[0],
                              salt_bridges=x.get('salt_bridges', [[]])[0], hbonds=x.get('hbonds', [[]])[0])
                cells.append(f"{x['resname'][:3].title()} {f(e, 1)}±{f(ese, 1)}（LJ {f(el, 1)}）；SASA {f(x['rel_sasa_mean'][0], 0)}%")
            else:
                cells.append('—')
        sysname = d.get('system')
        OUT['residue'][os.path.basename(p)] = dict(system=sysname, window=w, image_check=img, positions=rec)
        md.append(f"| {os.path.basename(p)} | {w[0]}–{w[1]}（{w[2]} 帧） | {img} | " + ' | '.join(cells) + ' |')
    md.append('')
    # spread across systems within each window (same residue at these positions in TUBA1A / 1B / 4A)
    by_w = {}
    for fn, r in OUT['residue'].items():
        if r['system'] in ('a1a', 'a1b', 'a4a'):
            by_w.setdefault(f"{r['window'][0]:.0f}-{r['window'][1]:.0f}", []).append(r)
    for wl, rs in sorted(by_w.items()):
        if len(rs) < 2:
            continue
        md.append(f"跨同种型范围（窗口 {wl} ns，{'/'.join(r['system'] for r in rs)}，单条轨迹）：\n")
        md.append('| 位点 | 残基 | 局部能量范围 kcal/mol | LJ 范围 | 相对 SASA 范围 % | 盐桥（占有率 ≥ 0.1） |')
        md.append('|---|---|---|---|---|---|')
        OUT.setdefault('residue_spread', {})[wl] = {}
        for k in POS:
            v = [r['positions'][k] for r in rs if k in r['positions']]
            if not v:
                continue
            E = [x['E_local'] for x in v]; LJ = [x['E_lj'] for x in v]; SA = [x['rel_sasa'] for x in v]
            sb = sorted({f'{b[0]}' for x in v for b in (x['salt_bridges'] or [])})
            OUT['residue_spread'][wl][k] = dict(E=[min(E), max(E)], LJ=[min(LJ), max(LJ)], SASA=[min(SA), max(SA)])
            md.append(f"| {k} | {v[0]['resname']} | {f(min(E), 1)} … {f(max(E), 1)} | {f(min(LJ), 1)} … {f(max(LJ), 1)} | "
                      f"{f(min(SA), 0)} … {f(max(SA), 0)} | {', '.join(sb) if sb else '—'} |")
        md.append('')
    # TUBA8 against TUBA1A / 1B / 4A on the same window: a difference counts only if TUBA8 lies outside the range of all
    # three by more than 3 combined block SEs of the nearest one (the criterion used since 2026-11-05 01:25)
    by_w8 = {}
    for fn, r in OUT['residue'].items():
        if r['system'] in ('a1a', 'a1b', 'a4a', 'a8'):
            by_w8.setdefault(f"{r['window'][0]:.0f}-{r['window'][1]:.0f}", {})[r['system']] = r
    for wl, rs in sorted(by_w8.items()):
        if 'a8' not in rs or len(rs) < 3:
            continue
        others = [rs[k] for k in ('a1a', 'a1b', 'a4a') if k in rs]
        md.append(f"TUBA8 与 {'/'.join(r['system'] for r in others)} 对比（窗口 {wl} ns；差值 = TUBA8 减去最接近的那个同种型；"
                  f"超出其余三者范围且大于 3 倍合并块误差才记为差异）：\n")
        md.append('| 位点 | TUBA8 | 其余 | 局部能量 TUBA8 / 其余范围 | 差值（倍合并 SE） | 其中 LJ TUBA8 / 其余范围 | 差值（倍合并 SE） | 相对 SASA TUBA8 / 其余范围 % | 差值（倍合并 SE） |')
        md.append('|---|---|---|---|---|---|---|---|---|')
        OUT.setdefault('residue_tuba8_vs_others', {})[wl] = {}
        for k in POS:
            x8 = rs['a8']['positions'].get(k)
            ov = [r['positions'][k] for r in others if k in r['positions']]
            if not x8 or not ov:
                continue
            def cmp(key, sekey):
                v8 = x8[key]; vals = [o[key] for o in ov]
                near = min(ov, key=lambda o: abs(o[key] - v8))
                se = math.hypot(x8.get(sekey) or 0, near.get(sekey) or 0) or 1e-9
                d = v8 - near[key]
                outside = v8 < min(vals) or v8 > max(vals)
                return v8, min(vals), max(vals), d, abs(d) / se, (outside and abs(d) > 3 * se)
            for o in [x8] + ov:
                o.setdefault('SASA_se', None)
            e = cmp('E_local', 'E_se')
            lj = cmp('E_lj', 'E_lj_se')
            sa = cmp('rel_sasa', 'SASA_se') if all(o.get('SASA_se') is not None for o in [x8] + ov) else None
            OUT['residue_tuba8_vs_others'][wl][k] = dict(E=e, LJ=lj, SASA=sa, res8=x8['resname'], res_other=ov[0]['resname'])
            md.append(f"| {k} | {x8['resname'].title()} | {ov[0]['resname'].title()} | {f(e[0], 1)} / {f(e[1], 1)} … {f(e[2], 1)} | "
                      f"{f(e[3], 1)}（{f(e[4], 1)}）{' **差异**' if e[5] else ''} | "
                      f"{f(lj[0], 1)} / {f(lj[1], 1)} … {f(lj[2], 1)} | {f(lj[3], 1)}（{f(lj[4], 1)}）{' **差异**' if lj[5] else ''} | "
                      f"{f(x8['rel_sasa'], 0)} / {f(min(o['rel_sasa'] for o in ov), 0)} … {f(max(o['rel_sasa'] for o in ov), 0)} | "
                      f"{(f(sa[3], 0) + '（' + f(sa[4], 1) + '）' + (' **差异**' if sa[5] else '')) if sa else '—'} |")
        md.append('')

# ---------------------------------------------------------------- 450 K unfolding
UNF = [(os.path.basename(p)[len('unfold_md450_'):-5], p) for p in sorted(glob.glob(os.path.join(I, 'unfold_md450_*.json')))]
UNF += [tuple(x.split('=', 1)) for x in a.unfold_extra]
for sid, p in UNF:
    d = load(p)
    OUT.setdefault('unfold', {})[sid] = d
    reps = sorted(d.get('replicas', {}))
    md.append(f'## 450 K 去折叠（a18，{sid}；Q_u = 该区段起始天然重原子接触的分数，t_half = 1 ns 滑动平均首次 < 0.5 的时间）\n')
    md.append('| 区段 | 残基 | 天然接触数 | ' + ' | '.join(f'{r}：t_half ns / 末 1 ns Q / 螺旋 起→末' for r in reps) + ' |')
    md.append('|---|---|---|' + '---|' * len(reps))
    for g, (lo, hi, n) in d['groups'].items():
        cells = []
        nn = None
        for r in reps:
            x = d['replicas'][r].get(g)
            if not x:
                cells.append('—'); continue
            nn = x['n_native_contacts']
            th = x['t_half_ns']
            T = d['replicas'][r]['ns']
            cells.append(f"{f(th, 1) if th is not None else f'> {T:.0f}（未达）'} / {f(x['Q_final_1ns'])} / "
                         f"{f(x.get('helix_start'))}→{f(x.get('helix_end'))}")
        md.append(f'| {g} | {lo}–{hi}（{n}） | {nn} | ' + ' | '.join(cells) + ' |')
    md.append('')
    for r in reps:
        R = d['replicas'][r]
        md.append(f"- {r}（{f(R['ns'], 1)} ns，{R['n_frames']} 帧）：去折叠顺序 {' → '.join(R.get('unfolding_order_units', []))}；"
                  f"GTP–蛋白接触原子数 起/末 {R.get('gtp_protein_contacts_first_last')}")
    md.append('')

if len(OUT.get('unfold', {})) >= 2:
    U = OUT['unfold']
    sids = list(U)
    groups = list(U[sids[0]]['groups'])
    md.append('450 K 各区段 t_half（ns；每个体系按副本列出，"—" 表示 20 ns 内未跌破 0.5）：\n')
    md.append('| 区段 | ' + ' | '.join(sids) + ' |')
    md.append('|---|' + '---|' * len(sids))
    for g in groups:
        cells = []
        for sid in sids:
            reps = sorted(U[sid]['replicas'])
            cells.append(' / '.join(f(U[sid]['replicas'][r].get(g, {}).get('t_half_ns'), 1) for r in reps))
        md.append(f'| {g} | ' + ' | '.join(cells) + ' |')
    md.append('')

# ---------------------------------------------------------------- REST2 progress
R2 = [(os.path.basename(p)[len('rest2_'):-5], p) for p in sorted(glob.glob(os.path.join(I, 'rest2_*.json')))]
R2 += [tuple(x.split('=', 1)) for x in a.rest2_extra]
if R2:
    md.append('## REST2（310 K 态；N/I/U 用同一校准阈值）\n')
    md.append('| 体系 | 每副本 ns | 往返 | g | 独立样本 | p(N) / p(I) / p(U) | 块 SD p(N) / p(U) | p(φ43 > 0) | 平均 Q(λ=1.2) / 环 RMSD Å |')
    md.append('|---|---|---|---|---|---|---|---|---|')
    OUT['rest2'] = {}
    for sid, p in R2:
        d = load(p)
        b = d.get('block_sd') or {}
        OUT['rest2'][sid] = {k: d.get(k) for k in ('p', 'block_sd', 'ns_per_replica', 'round_trips', 'g', 'n_samples',
                                                   'state0_direct', 'p_phi43_left', 'mean_Q_loop', 'mean_rmsd_loop_A')}
        md.append(f"| {sid} | {f(d.get('ns_per_replica'), 1)} | {d.get('round_trips')} | {f(d.get('g'), 0)} | {d.get('n_samples')} | "
                  f"{f(d['p']['N'], 3)} / {f(d['p']['I'], 3)} / {f(d['p']['U'], 3)} | {f(b.get('p_N'), 3)} / {f(b.get('p_U'), 3)} | "
                  f"{f(d.get('p_phi43_left'), 2)} | {f(d.get('mean_Q_loop'), 3)} / {f(d.get('mean_rmsd_loop_A'), 2)} |")
    md.append('')

# ---------------------------------------------------------------- tail binding (a22): ddG_bind and tauRAMD
BD = a.binding_dir or os.path.join(I, 'binding')
tb = os.path.join(BD, 'tail_rbfe.json')
if os.path.exists(tb):
    T = load(tb)
    OUT['tail_rbfe'] = T
    md.append('## 尾肽结合（a22；ΔΔG(Y→F) > 0 表示 Phe 尾结合更弱；kcal/mol）\n')
    md.append('| 量 | 值 | SD | 收敛门 | 计算方法 |')
    md.append('|---|---|---|---|---|')
    for k, v in T.items():
        md.append(f"| {v.get('quantity', k)} | {f(v.get('value'), 2 if 'Kd' not in k else 3)} | {f(v.get('sd'), 2)} | {v.get('gate', '')} | {v.get('method', '')} |")
    md.append('')
rr = os.path.join(BD, 'ramd_results.json')
if os.path.exists(rr):
    RM = load(rr)
    OUT['ramd'] = RM
    if RM:
        md.append('## τRAMD 停留时间（每条轨迹的解离时间取中位数，bootstrap 95 % 区间；未解离的按 20 ns 截断计，为下界）\n')
        md.append('| 体系 | τ ns | 95 % 区间 | 轨迹数 | 截断数 | koff 相对 TUBA1A 尾（EEY） |')
        md.append('|---|---|---|---|---|---|')
        for k, v in RM.items():
            md.append(f"| {k} | {f(v['tau_ns'], 2)} | {f(v['ci95'][0], 2)} – {f(v['ci95'][1], 2)} | {v['n']} | {v['censored']} | "
                      f"{f(v.get('koff_relative_to_Tuba1a_EEY'), 2)} |")
        md.append('')

# ---------------------------------------------------------------- PMF (a15) and MM/GBSA (a12b)
if a.pmf_dir:
    pf = sorted(glob.glob(os.path.join(a.pmf_dir, 'pmf', '*', 'pmf_results.json'))) + \
         sorted(glob.glob(os.path.join(a.pmf_dir, 'pmfdom', '*', 'pmf_results.json')))
    if pf:
        md.append('## 伞形采样 PMF（a15；dG0 用标准态 1 M；体相拟合残差 < 0.3 kcal/mol 才算达到体相；推出的 koff 来自一维 Smoluchowski，标为推导值）\n')
        md.append('| 体系 | dG0 | Kd M | 势垒 | r_min nm | 位点边界 nm | 体相达标（残差） | 最小相邻重叠 | 最小交换接受率 | 空箱 / 箱外样本 | koff s⁻¹（推导） |')
        md.append('|---|---|---|---|---|---|---|---|---|---|---|')
        OUT['pmf'] = {}
        for p_ in pf:
            d = load(p_); sid = '/'.join(p_.split(os.sep)[-3:-1])
            acc = d.get('acceptance') or [None]
            OUT['pmf'][sid] = {k: d.get(k) for k in ('dG0_kcal', 'Kd_M', 'barrier_kcal', 'r_min_nm', 'r_site_boundary_nm',
                               'gate_bulk_reached', 'bulk_fit_residual_kcal', 'min_adjacent_overlap', 'n_empty_bins',
                               'n_samples_outside_edges', 'koff_s')}
            ko = d.get('koff_s')
            md.append(f"| {sid} | {f(d.get('dG0_kcal'), 2)} | {d.get('Kd_M'):.2e} | {f(d.get('barrier_kcal'), 2)} | {f(d.get('r_min_nm'), 2)} | "
                      f"{f(d.get('r_site_boundary_nm'), 2)} | {'是' if d.get('gate_bulk_reached') else '否'}（{f(d.get('bulk_fit_residual_kcal'), 2)}） | "
                      f"{f(d.get('min_adjacent_overlap'), 3)} | {f(min(acc) if acc[0] is not None else None, 3)} | "
                      f"{d.get('n_empty_bins')} / {d.get('n_samples_outside_edges')} | {('%.2e' % ko) if isinstance(ko, (int, float)) and math.isfinite(ko) else '—'} |")
        md.append('')
if a.mmgbsa_dir:
    mf = sorted(glob.glob(os.path.join(a.mmgbsa_dir, '*.json')))
    if mf:
        md.append('## MM/GBSA（a12b；单轨迹、OBC2 mbondi2、0.15 M、不含熵；只作排序分数，不是亲和力；kcal/mol）\n')
        md.append('| 体系 | 副本均值 | 副本 SD | 副本数 | 各副本均值（帧数，重新成像帧数） |')
        md.append('|---|---|---|---|---|')
        OUT['mmgbsa'] = {}
        for p_ in mf:
            d = load(p_); sid = d.get('system', os.path.basename(p_)[:-5]); G = d.get('GB') or {}
            OUT['mmgbsa'][sid] = d
            reps = '；'.join(f"{k} {f(v['mean'], 1)}（{v['n_frames']}，{v.get('frames_ligand_reimaged', '—')}）" for k, v in sorted((d.get('replicas') or {}).items()))
            md.append(f"| {sid} | {f(G.get('mean'), 1)} | {f(G.get('sd_replicas'), 1)} | {G.get('n_replicas')} | {reps} |")
        md.append('')

# ---------------------------------------------------------------- US along phi43 (a17)
us = sorted(glob.glob(os.path.join(I, 'usphi_*.json')))
if us:
    md.append('## φ43 伞形采样（a17；W(φ) 由 MBAR 得到；ΔG = 左手（φ > 0）盆地减右手盆地，kcal/mol）\n')
    md.append('| 体系 | ΔG(左 − 右) | p(左手) | 经 φ = 0 的势垒 | 经 φ = 180° 的势垒 |')
    md.append('|---|---|---|---|---|')
    OUT['usphi'] = {}
    for p in us:
        d = load(p); sid = os.path.basename(p)[len('usphi_'):-5]
        OUT['usphi'][sid] = {k: d.get(k) for k in ('dG_left_minus_right_kcal', 'p_left', 'barrier_via_0_kcal', 'barrier_via_180_kcal')}
        md.append(f"| {sid} | {f(d.get('dG_left_minus_right_kcal'), 2)} | {f(d.get('p_left'), 3)} | {f(d.get('barrier_via_0_kcal'), 2)} | "
                  f"{f(d.get('barrier_via_180_kcal'), 2)} |")
    md.append('')

# ---------------------------------------------------------------- MSM (a11)
for spec in a.msm:
    lab, p = spec.split('=', 1)
    if not os.path.exists(p):
        continue
    d = load(p)
    OUT.setdefault('msm', {})[lab] = d
    md.append(f"## MSM（a11，{lab}；滞后 {d.get('lag_ns')} ns，TICA {d.get('tica_dim')} 维，{d.get('n_meta')} 个宏态）\n")
    its = d.get('implied_timescales_ns') or {}
    md.append('隐含时间尺度（最慢 3 个，ns）随滞后时间：' + '；'.join(f"{k} ns: {', '.join(f(x, 0) for x in v[:3])}" for k, v in its.items()) + '\n')
    md.append('| 体系 | CK 检验 | CK 最大偏差 | 宏态布居 | 最慢时间尺度 ns |')
    md.append('|---|---|---|---|---|')
    for sid, v in (d.get('systems') or {}).items():
        rec = v.get('failed_ck', v)
        md.append(f"| {sid} | {'通过' if rec.get('ck_pass') else ('未通过' if rec.get('ck_pass') is False else '没有跑成')} | "
                  f"{rec.get('ck_max_abs_dev') if isinstance(rec.get('ck_max_abs_dev'), str) else f(rec.get('ck_max_abs_dev'), 3)} | "
                  f"{', '.join(f(x, 3) for x in (rec.get('macro_populations') or []))} | "
                  f"{', '.join(f(x, 1) for x in (rec.get('slowest_timescales_ns') or []))} |")
    md.append('')

# ---------------------------------------------------------------- docking controls
if a.dock_dir:
    ds = sorted(glob.glob(os.path.join(a.dock_dir, '*.json')))
    if ds:
        md.append('## 对接对照（HADDOCK 3；参照 = AF3 复合物；分数是排序分数不是能量）\n')
        md.append('| 运行 | 模式 | 簇数 | 第 1 簇：n / HADDOCK 分数 / i-RMSD Å / Fnat / L-RMSD Å / DockQ |')
        md.append('|---|---|---|---|')
        OUT['dock'] = {}
        for p in ds:
            d = load(p)
            c = d['clusters'][0]
            OUT['dock'][os.path.basename(p)] = dict(mode=d['mode'], n_clusters=d['n_clusters'], top=c)
            md.append(f"| {os.path.basename(p)[:-5]} | {d['mode']} | {d['n_clusters']} | {c['n']} / {c['score']} / "
                      f"{c['irmsd']} / {c['fnat']} / {c['lrmsd']} / {c['dockq']} ± {c.get('dockq_std')} |")
        md.append('')

with open(a.out + '.md', 'w') as fh:
    fh.write('\n'.join(md) + '\n')
with open(a.out + '.json', 'w') as fh:
    json.dump(OUT, fh, indent=1, ensure_ascii=False, default=float)
print('wrote', a.out + '.md', a.out + '.json')
