"""Combine the tail alchemical legs into binding free-energy differences, and analyse the tauRAMD runs.

Tail (terminal residue Y<->F), from a09 leg summaries (fep_results.json):
  ddG_bind,TBCB(Y451F | Tuba1a tail) = dG(t_a1a_cplx Y451F) - dG(t_a1a_free Y451F)
  ddG_bind,TBCB(F449Y | Tuba8 tail)  = dG(t_a8_cplx F449Y)  - dG(t_a8_free F449Y)
  the same with TBCB K188A, and in the TBCE groove (e_*);
  coupling of the hydroxyl with Lys188: ddG_int = ddG(Y->F | K188A) - ddG(Y->F | WT)
  transfer selectivity TBCB -> TBCE: ddG_transfer(Y->F) = ddG_bind,TBCE(Y->F) - ddG_bind,TBCB(Y->F)
  Kd(F)/Kd(Y) = exp(ddG(Y->F)/kT)  (derived).
Positive ddG(Y->F) = the Phe tail binds more weakly.
tauRAMD: per system, tau = median dissociation time over trajectories, bootstrap (2000 resamples) 95 % interval;
relative koff = tau(reference)/tau(system).  Censored trajectories (no exit) are counted at the cap (lower bound).
Usage: a22_binding_analysis.py <results_dir>   (reads fep_results.json and ramd/*/s*/ramd_times.json)
"""
import os, sys, json, glob, math
import numpy as np
R = sys.argv[1]
kT = 0.0019872043 * 310.0
fr = json.load(open(os.path.join(R, 'fep_results.json')))['legs']
def leg(name):
    v = fr.get(name, {}).get('summary', {})
    return v.get('mean_dG_kcal'), v.get('sd_dG_kcal') or 0.0, v
out = {}
def diff(tag, a, b, sign=1, desc=''):
    (x, sx, vx), (y, sy, vy) = leg(a), leg(b)
    if x is None or y is None:
        out[tag] = dict(quantity=desc, value=None, sd=None, method='HREX/MBAR legs', prov='not_computed', gate=f'missing {a if x is None else b}')
        return None
    val = sign * (x - y); sd = math.hypot(sx, sy)
    gate = 'pass' if vx.get('gate_half_split') and vy.get('gate_half_split') else 'check'
    out[tag] = dict(quantity=desc, value=val, sd=sd, method=f'({a}) - ({b})', gate=gate)
    out[tag + '_Kd_ratio'] = dict(quantity=desc.replace('ddG', 'Kd ratio') + ' = exp(ddG/kT)', value=math.exp(val / kT),
                                  sd=None, method='derived from ddG', prov='derived', gate=gate)
    return val, sd
yf_wt = diff('TBCB_Y451F_Tuba1a_tail', 't_a1a_cplx_451PHE_folded', 't_a1a_free_451PHE_folded', 1, 'ddG_bind TBCB (Y451F, Tuba1a tail)')
fy_wt = diff('TBCB_F449Y_Tuba8_tail', 't_a8_cplx_449TYR_folded', 't_a8_free_449TYR_folded', 1, 'ddG_bind TBCB (F449Y, Tuba8 tail)')
yf_k = diff('TBCB-K188A_Y451F_Tuba1a_tail', 't_a1a_K188A_cplx_451PHE_folded', 't_a1a_free_451PHE_folded', 1, 'ddG_bind TBCB-K188A (Y451F, Tuba1a tail)')
fy_k = diff('TBCB-K188A_F449Y_Tuba8_tail', 't_a8_K188A_cplx_449TYR_folded', 't_a8_free_449TYR_folded', 1, 'ddG_bind TBCB-K188A (F449Y, Tuba8 tail)')
yf_e = diff('TBCE_Y451F_Tuba1a_tail', 'e_a1a_cplx_451PHE_folded', 't_a1a_free_451PHE_folded', 1, 'ddG_bind TBCE (Y451F, Tuba1a tail)')
fy_e = diff('TBCE_F449Y_Tuba8_tail', 'e_a8_cplx_449TYR_folded', 't_a8_free_449TYR_folded', 1, 'ddG_bind TBCE (F449Y, Tuba8 tail)')
def comb(tag, p, q, desc):
    if p is None or q is None:
        out[tag] = dict(quantity=desc, value=None, sd=None, method='', prov='not_computed', gate='missing leg'); return
    out[tag] = dict(quantity=desc, value=p[0] - q[0], sd=math.hypot(p[1], q[1]), method='difference of ddG', gate='')
comb('coupling_K188_hydroxyl_Tuba1a_tail', yf_k, yf_wt, 'ddG_int = ddG(Y->F | K188A) - ddG(Y->F | WT), Tuba1a tail')
if fy_k and fy_wt:
    out['coupling_K188_hydroxyl_Tuba8_tail'] = dict(quantity='ddG_int = ddG(Y->F | K188A) - ddG(Y->F | WT), Tuba8 tail (from F->Y legs)',
                                                    value=-(fy_k[0] - fy_wt[0]), sd=math.hypot(fy_k[1], fy_wt[1]), method='difference of ddG', gate='')
comb('transfer_selectivity_Y451F_Tuba1a', yf_e, yf_wt, 'ddG_transfer(Y->F) = ddG_bind,TBCE - ddG_bind,TBCB, Tuba1a tail')
if fy_e and fy_wt:
    out['transfer_selectivity_F449Y_Tuba8'] = dict(quantity='ddG_transfer(F->Y) = ddG_bind,TBCE - ddG_bind,TBCB, Tuba8 tail',
                                                   value=fy_e[0] - fy_wt[0], sd=math.hypot(fy_e[1], fy_wt[1]), method='difference of ddG', gate='')
json.dump(out, open(os.path.join(R, 'tail_rbfe.json'), 'w'), indent=1)
# ---- tauRAMD
ram = {}
for sd in sorted(glob.glob(os.path.join(R, 'ramd', '*'))):
    times = []; cens = 0
    for f in glob.glob(os.path.join(sd, 's*', 'ramd_times.json')):
        for t in json.load(open(f)):
            if t['t_ns'] is None:
                cens += 1; times.append(20.0)
            else:
                times.append(t['t_ns'])
    if len(times) < 5:
        continue
    x = np.array(times); rng = np.random.default_rng(1)
    boot = [np.median(rng.choice(x, len(x))) for _ in range(2000)]
    ram[os.path.basename(sd)] = dict(tau_ns=float(np.median(x)), ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
                                     n=len(x), censored=cens, gate='ok' if cens == 0 else f'{cens} censored (lower bound)')
ref = ram.get('t_a1a_cplx')
if ref:
    for k, v in ram.items():
        v['koff_relative_to_Tuba1a_EEY'] = ref['tau_ns'] / v['tau_ns']
json.dump(ram, open(os.path.join(R, 'ramd_results.json'), 'w'), indent=1)
print(json.dumps(out, indent=1)[:2000]); print(json.dumps(ram, indent=1)[:2000])
