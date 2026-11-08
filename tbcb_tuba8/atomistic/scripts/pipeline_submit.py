"""Wire the analysis stages to the compute jobs with Slurm dependencies (Sherlock).  Run on a login node; it only
submits / updates jobs (no polling).  Re-run it after resubmitting any compute job: pending stage jobs get their
dependency list rewritten to the current compute job ids, missing stage jobs are submitted.
Stage jobs run on the CPU partition 'normal' except mmgbsa (one GPU).  The report job depends on every compute job
and every stage job and mails the user (Slurm default address) on END or FAIL.
"""
import os, re, json, subprocess, sys
HOME = os.environ['HOME']; SCR = os.environ['SCRATCH']
J = f'{HOME}/tbcb_atomistic/jobs/pipeline'; os.makedirs(J, exist_ok=True)
LOG = f'{SCR}/tbcb_atomistic/logs'
STATE = f'{J}/pipeline.json'
GROUPS = {
    'fep': r'^fep_(a1a|a8)_114',
    'rest2': r'^rest2_',
    'md310': r'^md310_a',
    'md450': r'^md450_',
    'pmf': r'^pmf',
    'usphi': r'^usphi43_',
    'tail': r'^(fep_t_|fep_e_|ramd_)',
    'mmgbsa': r'^md310_(c_|x_|h_)',
}
q = subprocess.run(['squeue', '--me', '-h', '-o', '%i %j %T'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.split('\n')
jobs = [l.split() for l in q if l.strip()]
# only this campaign's compute jobs (2026-11-07: the report job had been made to wait for another project's jobs)
CAMPAIGN = re.compile('|'.join(GROUPS.values()))
mine = {jid: name for jid, name, st in jobs if CAMPAIGN.search(name)}
state = json.load(open(STATE)) if os.path.exists(STATE) else {}
stage_ids = {}
for st, pat in GROUPS.items():
    deps = sorted(j for j, n in mine.items() if re.search(pat, n))
    gpu = st == 'mmgbsa'
    hdr = ['sbatch', '--parsable', '-J', f'ana_{st}', '-o', f'{LOG}/ana_{st}_%j.out', '-t', '4:00:00', '--mail-type=FAIL']
    hdr += (['-p', 'gpu', '-G', '1', '-c', '4', '--mem=24G'] if gpu else ['-p', 'normal', '-c', '4', '--mem=16G'])
    dep = 'afterany:' + ':'.join(deps) if deps else ''
    old = state.get(st)
    alive = old and old in [j for j, n, s in jobs]
    if alive:
        if dep:
            subprocess.run(['scontrol', 'update', f'JobId={old}', f'Dependency={dep}'], check=False)
        stage_ids[st] = old
        print(f'{st}: updated job {old} -> {len(deps)} compute deps')
    else:
        # a stage that already ran (results/stage_<st>.json) and has no compute job left to wait for is not run again
        rec = f'{SCR}/tbcb_atomistic/results/stage_{st}.json'
        if not deps and os.path.exists(rec) and not json.load(open(rec)).get('failed'):
            print(f'{st}: ran earlier without failures and nothing new to wait for; not resubmitted')
            continue
        cmd = hdr + ([f'--dependency={dep}'] if dep else []) + [f'{HOME}/tbcb_atomistic/scripts/pipeline_stage.sh', st]
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        jid = r.stdout.strip().split(';')[0]
        if not jid:
            print(st, 'submit failed:', r.stderr.strip()); continue
        stage_ids[st] = jid
        print(f'{st}: submitted {jid} after {len(deps)} compute jobs')
state.update(stage_ids)
alldeps = sorted(set(mine) | set(stage_ids.values()))
rep = state.get('report')
alive = rep and rep in [j for j, n, s in jobs]
dep = 'afterany:' + ':'.join(alldeps)
if alive:
    subprocess.run(['scontrol', 'update', f'JobId={rep}', f'Dependency={dep}'], check=False)
    print('report: updated', rep, len(alldeps), 'deps')
else:
    r = subprocess.run(['sbatch', '--parsable', '-J', 'report_tbcb', '-o', f'{LOG}/report_%j.out', '-p', 'normal', '-c', '2',
                        '--mem=8G', '-t', '2:00:00', '--mail-type=END,FAIL', f'--dependency={dep}',
                        f'{HOME}/tbcb_atomistic/scripts/pipeline_stage.sh', 'report'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    state['report'] = r.stdout.strip().split(';')[0]
    print('report: submitted', state['report'], 'after', len(alldeps), 'jobs', r.stderr.strip())
json.dump(state, open(STATE, 'w'), indent=1)
print(json.dumps(state))
