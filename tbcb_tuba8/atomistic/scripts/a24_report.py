"""Write the Chinese results report (Markdown) from the master table (a21_compile.py output) plus the force table.
Every number carries its provenance class (computed / derived / score / not_computed) and its gate status.
Usage: a24_report.py <results_dir> <out.md>"""
import os, sys, json, csv, math
R, out = sys.argv[1], sys.argv[2]
rows = json.load(open(os.path.join(R, 'master.json')))
PROV = {'computed': '实算物理量', 'derived': '由实算量经物理模型推导', 'score': '排序分数（非自由能，不可换算Kd）',
        'not_computed': '未计算'}
BLOCK = [('structures', '一、起始结构的能量最小化与残余原子受力'), ('folding', '二、残基替换的折叠自由能 ΔΔG_fold'),
         ('fep_leg', '二(附)、各炼金臂的自由能'), ('intermediates', '三、折叠中间态：布居、自由能面与能垒'),
         ('kinetics', '四、动力学：MSM、MFPT、解离速率与停留时间'), ('binding', '五、TBCB/TBCE 结合：ΔΔG、ΔG°、Kd 与打分'),
         ('not_computed', '六、要求了但没有算出的量及原因')]
def fmt(v):
    if v is None or v == '':
        return '—'
    try:
        x = float(v)
    except Exception:
        return str(v)
    if x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5):
        return f'{x:.3e}'
    return f'{x:.3f}'
L = ['# TBCB–TUBA8 原子级中间态物理计算：结果总表', '',
     '本表只收录原子级模拟直接得到的量，以及由这些量经明确物理模型推导的量；不含任何动力学网络模型拟合。',
     '每一行标注：实算物理量 / 推导量 / 排序分数 / 未计算，以及该数是否通过协议规定的验证门。', '']
for b, title in BLOCK:
    sel = [r for r in rows if r['block'] == b]
    if not sel:
        continue
    L += [f'## {title}', '', '| 体系 | 物理量 | 数值 | 不确定度 | 单位 | 方法 | 类别 | 验证门 |', '|---|---|---|---|---|---|---|---|']
    for r in sel:
        unc = r['uncertainty']
        if isinstance(unc, list):
            unc = '–'.join(fmt(u) for u in unc)
        else:
            unc = fmt(unc)
        L.append(f"| {r['system']} | {r['quantity']} | {fmt(r['value'])} | {unc} | {r['unit']} | {r['method']} | "
                 f"{PROV.get(r['provenance'], r['provenance'])} | {r['gate']} |")
    L.append('')
open(out, 'w').write('\n'.join(L) + '\n')
print('wrote', out, len(rows), 'rows')
