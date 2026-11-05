"""Write one sbatch file per job for Sherlock (gpu partition, 1 GPU, 47 h) and a submission list.
Each job runs one or more commands sequentially; every command resumes from its own checkpoint, so a job that
is cut short is continued by resubmitting the same file (no duplicates are ever run at once).
Usage: make_jobs.py <wave> ; writes $HOME/tbcb_atomistic/jobs/<wave>/*.sbatch and submit_<wave>.txt"""
import os, sys
HOME = os.environ.get('HOME'); SCR = os.environ.get('SCRATCH')
B = f'{SCR}/tbcb_atomistic/build'; R = f'{SCR}/tbcb_atomistic/runs'
HDR = """#!/bin/bash
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t {hours}:00:00
#SBATCH -J {name}
#SBATCH -o {logs}/%x_%j.out
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
date
"""
WAVES = {}
def fep(sid, res, to, leg, rep, ns=5.0):
    return (f'fep_{sid}_{res}{to}_{leg}_r{rep}',
            f'scripts/a05_fep.py {B}/{sid} {R}/fep/{sid}_{res}{to}_{leg} {rep} --res {res} --to {to} --leg {leg} --ns {ns}')
def rest2(sid, rep=0, ns=30.0, n=12, hot='30-50'):
    return (f'rest2_{sid}_r{rep}', f'scripts/a04_rest2.py {B}/{sid} {R}/rest2/{sid} {rep} --hot {hot} --n {n} --ns {ns}')
def md(sid, rep, ns=200.0, T=310):
    extra = '' if T == 310 else f' --T {T} --nvt --save_ps 10'
    return (f'md{T}_{sid}_r{rep}', f'scripts/a03_md.py {B}/{sid} {R}/md{T}/{sid} {rep} {ns}{extra}')

w1 = []
for rep in range(3):
    w1.append([fep('a1a', 114, 'SER', 'folded', rep)])
    w1.append([fep('a8', 114, 'ILE', 'folded', rep)])
w1.append([fep('a1a', 114, 'SER', 'unfolded', r) for r in range(3)] + [fep('a8', 114, 'ILE', 'unfolded', r) for r in range(3)])
for sid in ('a1a', 'a8', 'a8_loop1A', 'a1a_loop8'):
    w1.append([rest2(sid)])
for sid in ('a1a', 'a8', 'a1b', 'a4a'):
    w1.append([md(sid, 0)])
w1.append([md(s, r, ns=20.0, T=450) for s in ('a1a', 'a8') for r in range(3)])
WAVES['w1'] = w1

w2 = []
for sid in ('a1a', 'a8', 'a1b', 'a4a'):
    for rep in (1, 2):
        w2.append([md(sid, rep)])
for sid in ('a1a_I114S', 'a8_S114I', 'a8_loop1A', 'a1a_loop8', 'a8_A40K', 'a1a_K40A', 'a8_I43G', 'a1a_G43I'):
    w2.append([md(sid, 0)])
for sid in ('a1b', 'a4a', 'a8_I43G', 'a1a_G43I', 'a8_A40K', 'a1a_K40A'):
    w2.append([rest2(sid)])
for res, to in ((110, 'VAL'), (294, 'SER'), (270, 'VAL')):
    for rep in range(2):
        w2.append([fep('a1a', res, to, 'folded', rep)])
    w2.append([fep('a1a', res, to, 'unfolded', r) for r in range(2)])
for rep in range(3):
    w2.append([fep('a1a_loop8', 114, 'SER', 'folded', rep)])
    w2.append([fep('a8_loop1A', 114, 'ILE', 'folded', rep)])
w2.append([md(s, r, ns=20.0, T=450) for s in ('a1b', 'a4a', 'a1a_I114S', 'a8_S114I') for r in range(3)])
WAVES['w2'] = w2

w3 = []
LIB = f'{B}/a1a/minimized.cif.gz'
def tailfep(sid, res, to, rep, restr):
    extra = f' --donor_lib {LIB}' + (f' --restrain "{restr}"' if restr else '')
    return (f'fep_{sid}_{res}{to}_r{rep}',
            f'scripts/a05_fep.py {B}/{sid} {R}/fep/{sid}_{res}{to}_folded {rep} --res {res} --to {to} --leg folded --ns 5.0{extra}')
for sid, res, to in (('t_a1a_cplx', 451, 'PHE'), ('t_a8_cplx', 449, 'TYR'), ('t_a1a_K188A_cplx', 451, 'PHE'), ('t_a8_K188A_cplx', 449, 'TYR')):
    for rep in range(3):
        w3.append([tailfep(sid, res, to, rep, f'A:187,188,203,204,223|B:{res}')])
w3.append([tailfep('t_a1a_free', 451, 'PHE', r, None) for r in range(3)])
w3.append([tailfep('t_a8_free', 449, 'TYR', r, None) for r in range(3)])
for sid, nrep in (('c_tbcb_a1a_m0', 3), ('c_tbcb_a8_m0', 3), ('c_tbce_a1a_m0', 2), ('c_tbce_a8_m0', 2),
                  ('x_9M1L_tbce_a1a', 2), ('x_9M1L_tbce_a8', 2)):
    for rep in range(nrep):
        w3.append([md(sid, rep, ns=50.0)])
WAVES['w3'] = w3

w3b = []
for sid, res, to in (('e_a1a_cplx', 451, 'PHE'), ('e_a8_cplx', 449, 'TYR')):
    for rep in range(3):
        w3b.append([tailfep(sid, res, to, rep, f'A:31,32,48,49,68,69|B:{res}')])
def ramd(sid, start, ntraj=10):
    return (f'ramd_{sid}_s{start}', f'scripts/a06_ramd.py {B}/{sid} {R}/ramd/{sid} {start} --ntraj {ntraj}')
RAMD = ['t_a1a_cplx', 't_a1aYF_cplx', 't_a8_cplx', 't_a8FY_cplx', 't_a1a_K188A_cplx', 't_a1aYF_K188A_cplx',
        't_a8_K188A_cplx', 't_a8FY_K188A_cplx', 'e_a1a_cplx', 'e_a8_cplx']
for sid in RAMD:
    w3b.append([ramd(sid, st) for st in range(3)])
WAVES['w3b'] = w3b

w3c = []
for sid in ('x_9M1L_tbce_a1a', 'x_9M1L_tbce_a8', 'h_9M1M_Tuba1a', 'h_9M1M_Tuba8', 'h_9M1M_Tuba1b', 'h_9M1M_Tuba4a'):
    for rep in range(2):
        w3c.append([md(sid, rep, ns=50.0)])
WAVES['w3c'] = w3c

w4 = []
for sid, groove, term in (('t_a1a_cplx_big', 'A:187,188,203,204,223', 'B:451'), ('t_a8_cplx_big', 'A:187,188,203,204,223', 'B:449'),
                          ('e_a1a_cplx_big', 'A:31,32,48,49,68,69', 'B:451'), ('e_a8_cplx_big', 'A:31,32,48,49,68,69', 'B:449')):
    w4.append([(f'pmf_{sid}_r0', f'scripts/a15_pmf_tail.py run {B}/{sid} {R}/pmf/{sid} 0 --groove {groove} --term {term} --ns 8')])
for sid in ('a1a', 'a8', 'a8_I43G', 'a1a_G43I', 'a8_loop1A', 'a1a_loop8'):
    w4.append([(f'usphi43_{sid}_r0', f'scripts/a17_us_phi.py run {B}/{sid} {R}/usphi43/{sid} 0 --res 43 --ns 4')])
WAVES['w4'] = w4

w5 = []
for sid in ('a1a_dom', 'a8_dom'):
    w5.append([(f'pmfdom_{sid}_r0', f'scripts/a15_pmf_tail.py run {B}/{sid} {R}/pmfdom/{sid} 0 --groove "*:1-205" --term "*:206-381" '
                                    f'--g1names CA --g2names CA --span 2.5 --ns 6')])
WAVES['w5'] = w5

wave = sys.argv[1]
J = f'{HOME}/tbcb_atomistic/jobs/{wave}'; L = f'{SCR}/tbcb_atomistic/logs'
os.makedirs(J, exist_ok=True); os.makedirs(L, exist_ok=True)
lst = []
for cmds in WAVES[wave]:
    name = cmds[0][0] if len(cmds) == 1 else cmds[0][0] + f'_x{len(cmds)}'
    hours = 47
    budget = 47 * 60 - 25
    body = HDR.format(hours=hours, name=name, logs=L)
    if len(cmds) == 1:
        body += f'python3 -W ignore {cmds[0][1]} --budget {budget}\n'
    else:
        body += 'T0=$(date +%s)\n'
        for n, c in cmds:
            body += (f'LEFT=$(( {budget} - ($(date +%s) - T0) / 60 ))\n'
                     f'[ $LEFT -gt 30 ] && python3 -W ignore {c} --budget $LEFT\n')
    body += 'date\n'
    fn = os.path.join(J, f'{name}.sbatch')
    open(fn, 'w').write(body)
    lst.append(fn)
open(os.path.join(J, f'submit_{wave}.txt'), 'w').write('\n'.join(lst) + '\n')
print(f'{len(lst)} jobs in {J}')
