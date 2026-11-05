"""Per-job status for this campaign (Python 3.6 compatible; read-only: one squeue + one sacct call, then logs).
Prints one JSON object per job: id, name, state, reason, elapsed, gpu, progress.
Usage: status_report.py <logdir> [--marlowe]"""
import os, re, sys, json, glob, subprocess
logdir = sys.argv[1]
def run(cmd):
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout
q = run(['squeue', '--me', '-h', '-o', '%i|%j|%T|%r|%M'] if '--marlowe' not in sys.argv else
        ['squeue', '-u', os.environ.get('USER', ''), '-h', '-o', '%i|%j|%T|%r|%M'])
rows = {}
for l in q.strip().split('\n'):
    if not l.strip():
        continue
    i, n, s, r, m = l.split('|')
    rows[i] = dict(id=i, name=n, state=s, reason=r, elapsed=m)
acct = run(['sacct', '-S', '2026-11-04T15:00', '-X', '-n', '-P', '-o', 'JobID,JobName,State,Elapsed'])
for l in acct.strip().split('\n'):
    if not l.strip():
        continue
    i, n, s, e = l.split('|')[:4]
    if i not in rows and not n.startswith('bash') and n not in ('zl_gputest',):
        rows[i] = dict(id=i, name=n, state=s.split()[0], reason='', elapsed=e)
pat = re.compile(r'^(it \d+/\d+.*|[\d.]+ ns +\d+ ns/day|segment \d+ done.*|traj \d+:.*|pull: window.*|sequential warm.*|'
                 r'iteration \d+/\d+.*|production on.*|equilibrating.*|== .*|HADDOCK_EXIT=.*|\{"stage".*\})$')
for i, r in rows.items():
    fs = glob.glob(os.path.join(logdir, '*_%s.out' % i))
    prog, gpu, err = '', '', ''
    if fs:
        try:
            txt = open(fs[0], errors='ignore').read().split('\n')
        except Exception:
            txt = []
        if txt and ('NVIDIA' in txt[0] or 'GeForce' in txt[0] or 'H100' in txt[0]):
            gpu = txt[0].split(',')[0]
        for l in txt:
            if pat.match(l.strip()):
                prog = l.strip()[:120]
            if 'Traceback' in l or 'Error' in l or 'NaN' in l:
                err = l.strip()[:120]
        if err and 'stage' not in prog:
            prog = prog + ' | ERR: ' + err
    r['progress'] = prog; r['gpu'] = gpu
    print(json.dumps(r))
