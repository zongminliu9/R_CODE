"""Collect every analysis JSON under results/ into one master table with provenance.

Provenance classes (the user's three categories):
  computed     - a physical quantity computed directly from the atomistic simulations here (FEP/MBAR free energies,
                 PMF-derived dG0/Kd, REST2/MBAR populations, MSM populations/MFPT, tauRAMD residence-time ratios,
                 minimisation energies/forces, geometric/interaction-energy statistics).
  derived      - obtained from computed quantities through an explicit physical model with stated assumptions
                 (e.g. koff/kon from a 1-D Smoluchowski treatment of a computed PMF; Kd ratios from ddG).
  score        - a ranking score that is not a free energy (MM/GBSA, MM/PBSA, HADDOCK score).  Never a Kd.
  not_computed - requested but not obtained, with the reason.
Each row: block, system, quantity, value, uncertainty, unit, method, provenance, gate, source.
Usage: a21_compile.py <results_dir> <out_prefix>
"""
import os, sys, json, glob, csv, math
R, out = sys.argv[1], sys.argv[2]
rows = []
def add(block, system, quantity, value, unc=None, unit='', method='', prov='computed', gate='', source=''):
    rows.append(dict(block=block, system=system, quantity=quantity, value=value, uncertainty=unc, unit=unit,
                     method=method, provenance=prov, gate=gate, source=os.path.relpath(source, R) if source else ''))
kT = 0.0019872043 * 310.0
# ---- folding ddG
f = os.path.join(R, 'fep_results.json')
if os.path.exists(f):
    d = json.load(open(f))
    for k, v in d.get('ddG_fold', {}).items():
        g = 'pass' if v['folded'].get('gate_half_split') and v['unfolded'].get('gate_half_split') else 'check'
        add('folding', k, 'ddG_fold', v['ddG_fold_kcal'], v['sd_kcal'], 'kcal/mol', 'dual-topology FEP, HREX, MBAR; folded - pentapeptide', gate=g, source=f)
        add('folding', k, 'folded/unfolded population ratio change at 310 K', math.exp(v['ddG_fold_kcal'] / kT), None, 'fold', 'exp(ddG/RT)', 'derived', source=f)
    for k, v in d.get('legs', {}).items():
        s = v['summary']
        add('fep_leg', k, 'dG_leg', s.get('mean_dG_kcal'), s.get('sd_dG_kcal'), 'kcal/mol', 'HREX/MBAR', gate=f"{s.get('n_passing')}/{s.get('n_replicates')} replicas pass", source=f)
# ---- binding RBFE (tail), double-mutant cycle, transfer selectivity
f = os.path.join(R, 'tail_rbfe.json')
if os.path.exists(f):
    d = json.load(open(f))
    for k, v in d.items():
        add('binding', k, v['quantity'], v['value'], v.get('sd'), 'kcal/mol', v['method'], v.get('prov', 'computed'), v.get('gate', ''), f)
# ---- PMF absolute binding
for f in sorted(glob.glob(os.path.join(R, 'pmf', '*', 'pmf_results.json'))):
    d = json.load(open(f)); sid = f.split(os.sep)[-2]
    gate = 'pass' if d.get('gate_bulk_reached') and d.get('min_adjacent_overlap', 0) >= 0.03 else 'fail'
    add('binding', sid, 'dG0_bind', d['dG0_kcal'], None, 'kcal/mol', 'US/HREX PMF + spherical standard state', gate=gate, source=f)
    add('binding', sid, 'Kd', d['Kd_M'], None, 'M', 'from dG0', 'derived', gate=gate, source=f)
    add('binding', sid, 'unbinding barrier (PMF)', d['barrier_kcal'], None, 'kcal/mol', 'PMF', gate=gate, source=f)
    for q in ('koff_s', 'kon_M_s', 'residence_time_s'):
        if q in d:
            add('kinetics', sid, q, d[q], None, '', '1-D Smoluchowski on PMF with D(r)', 'derived', gate=gate, source=f)
# ---- tauRAMD
f = os.path.join(R, 'ramd_results.json')
if os.path.exists(f):
    d = json.load(open(f))
    for k, v in d.items():
        add('kinetics', k, 'tauRAMD residence time (relative)', v['tau_ns'], v.get('ci95'), 'ns (RAMD time)', 'tauRAMD', gate=v.get('gate', ''), source=f)
# ---- REST2
for f in sorted(glob.glob(os.path.join(R, 'rest2', '*', 'rest2_results.json'))):
    d = json.load(open(f)); sid = f.split(os.sep)[-2]; sd = d.get('block_sd', {})
    for s in ('N', 'I', 'U'):
        add('intermediates', sid, f'p_loop_{s} (310 K)', d['p'][s], sd.get(f'p_{s}'), '', 'REST2 + MBAR reweighting', source=f)
    for s in ('I', 'U'):
        add('intermediates', sid, f'dG_loop_{s}-N', d['dG_kcal'][s], sd.get(f'dG_{s}'), 'kcal/mol', 'REST2 + MBAR', source=f)
    add('intermediates', sid, 'p(phi43>0)', d['p_phi43_left'], sd.get('p_phi43_left'), '', 'REST2 + MBAR', source=f)
    if d.get('barrier_N_to_U_kcal') is not None:
        add('intermediates', sid, 'loop N->U barrier along Q_loop', d['barrier_N_to_U_kcal'], None, 'kcal/mol', 'REST2 PMF', source=f)
# ---- phi43 umbrella
for f in sorted(glob.glob(os.path.join(R, 'usphi43', '*', 'usphi_results.json'))):
    d = json.load(open(f)); sid = f.split(os.sep)[-2]
    add('intermediates', sid, 'dG(phi43 left - right)', d['dG_left_minus_right_kcal'], None, 'kcal/mol', 'US/HREX PMF', source=f)
    add('intermediates', sid, 'phi43 barrier via 0', d['barrier_via_0_kcal'], None, 'kcal/mol', 'US/HREX PMF', source=f)
    add('intermediates', sid, 'phi43 barrier via 180', d['barrier_via_180_kcal'], None, 'kcal/mol', 'US/HREX PMF', source=f)
# ---- MSM
f = os.path.join(R, 'msm', 'msm_results.json')
if os.path.exists(f):
    d = json.load(open(f))
    for sid, v in d['systems'].items():
        if 'failed_ck' in v:
            add('kinetics', sid, 'MSM', None, None, '', 'MSM', 'not_computed', 'CK test failed: not reported', f); continue
        for i, p in enumerate(v['macro_populations']):
            add('kinetics', sid, f'MSM macrostate {i} population', p, None, '', f"MSM lag {d['lag_ns']} ns", gate='CK pass', source=f)
        for k, t in v['mfpt_ns'].items():
            add('kinetics', sid, f'MFPT {k}', t, None, 'ns', 'MSM', gate='CK pass', source=f)
# ---- domain-separation PMF
for f in sorted(glob.glob(os.path.join(R, 'pmfdom', '*', 'pmf_results.json'))):
    d = json.load(open(f)); sid = f.split(os.sep)[-2]
    W = [x for x in d['W_kcal'] if x is not None]
    add('intermediates', sid, 'N/intermediate-domain separation PMF, W(max) - W(min)', max(W) - min(W) if W else None, None, 'kcal/mol', 'US/HREX PMF', source=f)
# ---- MM/GBSA, MM/PBSA
for f in sorted(glob.glob(os.path.join(R, 'mmpbsa', '*.json'))):
    d = json.load(open(f)); sid = os.path.basename(f)[:-5]
    for k in ('GB', 'PB'):
        if k in d:
            add('binding', sid, f'MM/{k}SA dG (score)', d[k]['mean'], d[k].get('sd_replicas'), 'kcal/mol', f'MMPBSA.py {k}', 'score', 'replica spread', f)
# ---- minimisation forces
f = os.path.join(R, 'minimization_forces.json')
if os.path.exists(f):
    for r in json.load(open(f)):
        for q in ('init_PE_kcal', 'init_Fmax_kcal_A', 'init_RMS_solute_kcal_A', 'final_PE_kcal', 'final_Fmax_kcal_A',
                  'final_RMS_solute_kcal_A', 'vac_PE_kcal', 'vac_Fmax_kcal_A', 'vac_RMS_kcal_A', 'if_CA_rmsd_A'):
            if r.get(q) not in (None, ''):
                add('structures', r['system'], q, float(r[q]), None, '', 'OpenMM ff14SB/Meagher GTP/IOD Mg', source=f)
# ---- explicitly not computed
for sid, q, why in (
        ('TBCB:alpha full length', 'absolute dG_bind / Kd', '~5000 A^2 interface; no converged absolute method within reach; tail-site dG0 computed instead'),
        ('TBCB->TBCE transfer', 'reaction path, barrier, rate', 'no atomic structure of a TBCB:TBCE:alpha complex (EMD-2447 is 21 A); only end-state groove thermodynamics computed'),
        ('K40A / A40K', 'alchemical ddG_fold', 'charge-changing transformation; finite-size PME artefact larger than the signal; studied by MD/REST2 populations instead'),
        ('G43I / I43G', 'alchemical ddG_fold', 'glycine backbone parameters differ; studied by the phi43 umbrella PMF instead')):
    add('not_computed', sid, q, None, None, '', '', 'not_computed', why)
keys = ['block', 'system', 'quantity', 'value', 'uncertainty', 'unit', 'method', 'provenance', 'gate', 'source']
with open(out + '.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(rows, open(out + '.json', 'w'), indent=1, default=str)
print(len(rows), 'rows ->', out + '.csv')
