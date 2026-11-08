"""Marlowe (preempt, billed $0.25/GPU-h) job chains for wave 2.  Each simulation is a chain of 6-hour jobs linked
by afterany dependencies; every job resumes from the simulation's own checkpoint and exits at once if it is done.
OpenCL platform (the conda-forge CUDA build fails on these nodes with CUDA_ERROR_UNSUPPORTED_PTX_VERSION).
Writes jobs/<wave>/<name>.sbatch and chains_<wave>.txt (name, n_links)."""
import os, sys
ROOT = f"/scratch/m000251/{os.environ['USER']}/tbcb_mech/atomistic"
B = f'{ROOT}/build'; R = f'{ROOT}/runs'
PY = f"/scratch/m000251/{os.environ['USER']}/mm/envs/amb/bin/python"
HDR = """#!/bin/bash
#SBATCH -p preempt
#SBATCH -A marlowe-m000251
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 6:00:00
#SBATCH --exclude=n18,n28
#SBATCH -J {name}
#SBATCH -o {ROOT}/logs/%x_%j.out
nvidia-smi --query-gpu=name --format=csv,noheader
cd {ROOT}
date
"""
def fep(sid, res, to, leg, rep, ns=5.0):
    return (f'fep_{sid}_{res}{to}_{leg}_r{rep}', f'scripts/a05_fep.py {B}/{sid} {R}/fep/{sid}_{res}{to}_{leg} {rep} --res {res} --to {to} --leg {leg} --ns {ns}', 2)
def rest2(sid, rep=0, ns=30.0, n=12):
    return (f'rest2_{sid}_r{rep}', f'scripts/a04_rest2.py {B}/{sid} {R}/rest2/{sid} {rep} --hot 30-50 --n {n} --ns {ns}', 3)
def md(sid, rep, ns=200.0, T=310):
    extra = '' if T == 310 else f' --T {T} --nvt --save_ps 10'
    return (f'md{T}_{sid}_r{rep}', f'scripts/a03_md.py {B}/{sid} {R}/md{T}/{sid} {rep} {ns}{extra}', 2 if T == 310 else 1)
W = {}
w2 = []
for sid in ('a1a', 'a8'):
    for rep in (1, 2):
        w2.append([md(sid, rep)])
for sid in ('a1a_I114S', 'a8_S114I', 'a8_loop1A', 'a1a_loop8', 'a8_A40K', 'a1a_K40A'):
    w2.append([md(sid, 0)])
for sid in ('a1b', 'a4a', 'a8_I43G', 'a1a_G43I'):
    w2.append([rest2(sid)])
for res, to in ((110, 'VAL'), (294, 'SER'), (270, 'VAL')):
    for rep in range(2):
        w2.append([fep('a1a', res, to, 'folded', rep)])
    w2.append([fep('a1a', res, to, 'unfolded', r) for r in range(2)])
for rep in range(3):
    w2.append([fep('a1a_loop8', 114, 'SER', 'folded', rep)])
    w2.append([fep('a8_loop1A', 114, 'ILE', 'folded', rep)])
w2.append([md(s, r, ns=20.0, T=450) for s in ('a1b', 'a4a') for r in range(3)])
w2.append([md(s, r, ns=20.0, T=450) for s in ('a1a_I114S', 'a8_S114I') for r in range(3)])
W['w2'] = w2
wave = sys.argv[1]
J = f'{ROOT}/jobs/{wave}'; os.makedirs(J, exist_ok=True); os.makedirs(f'{ROOT}/logs', exist_ok=True)
rows = []
for cmds in W[wave]:
    name = cmds[0][0] if len(cmds) == 1 else cmds[0][0] + f'_x{len(cmds)}'
    links = max(c[2] for c in cmds) if len(cmds) == 1 else 2
    budget = 6 * 60 - 20
    body = HDR.format(name=name, ROOT=ROOT) + 'T0=$(date +%s)\n'
    for n, c, _ in cmds:
        body += (f'LEFT=$(( {budget} - ($(date +%s) - T0) / 60 ))\n'
                 f'[ $LEFT -gt 20 ] && {PY} -W ignore {c} --budget $LEFT --platform OpenCL\n')
    body += 'date\n'
    open(os.path.join(J, f'{name}.sbatch'), 'w').write(body)
    rows.append(f'{name}\t{links}')
open(os.path.join(J, f'chains_{wave}.txt'), 'w').write('\n'.join(rows) + '\n')
print(len(rows), 'chains,', sum(int(r.split()[1]) for r in rows), 'job links')
