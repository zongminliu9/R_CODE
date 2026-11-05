"""Table of energy minimisation and residual atomic forces for every starting structure in build/.
Columns: as-built (solvated, flexible system) potential energy, max force (atom), RMS force over all atoms and
over solute; the same after restrained and after free minimisation; backbone RMSD start->minimised; Mg2+ inner
shell after minimisation; interface heavy-atom contacts (4.5 A) and interface residues before/after and interface CA
RMSD for complexes; and for the AF3 / experimental complexes the vacuum as-built record (protein + GTP.Mg, no
solvent, no cutoff, hydrogens placed by OpenMM, before any relaxation).
Usage: a20_force_table.py <build_root> <out_prefix>"""
import os, sys, json, glob, csv
root, out = sys.argv[1], sys.argv[2]
rows = []
for d in sorted(glob.glob(os.path.join(root, '*'))):
    sid = os.path.basename(d)
    r = {'system': sid}
    b = os.path.join(d, 'build.json'); m = os.path.join(d, 'minimize.json'); v = os.path.join(d, 'vacuum_forces.json')
    if os.path.exists(b):
        bj = json.load(open(b))
        r['n_atoms'] = bj.get('n_atoms'); r['changes'] = ';'.join(bj.get('changes', [])) if isinstance(bj.get('changes'), list) else ''
        g = bj.get('grafts') or {}
        for k, gi in g.items():
            if 'min_dist_graft_to_partner_A' in gi:
                r[f'GTPgraft_{k}_to_partner_A'] = gi['min_dist_graft_to_partner_A']
        if 'built_loop_relaxation' in bj:
            r['rebuilt_gap_residues'] = bj['built_loop_relaxation']['n_residues']
        if bj.get('graft_site_ca_rmsd_A') is not None:
            r['graft_site_CA_rmsd_A'] = bj['graft_site_ca_rmsd_A']
    if os.path.exists(v):
        vj = json.load(open(v))
        r.update(vac_PE_kcal=round(vj['potential_kcalmol'], 1), vac_Fmax_kcal_A=round(vj['max_force_kcalmol_A'], 2),
                 vac_Fmax_atom=vj['max_force_atom'], vac_RMS_kcal_A=round(vj['rms_force_kcalmol_A'], 3),
                 vac_n_atoms_F_gt_100=vj['n_atoms_force_gt_100'])
    if os.path.exists(m):
        mj = json.load(open(m))
        for tag, key in (('init', 'initial'), ('stage1', 'after_stage1'), ('final', 'final')):
            s = mj[key]
            r[f'{tag}_PE_kcal'] = round(s['potential_kcalmol'], 1); r[f'{tag}_Fmax_kcal_A'] = round(s['max_force_kcalmol_A'], 2)
            r[f'{tag}_Fmax_atom'] = s['max_force_atom']; r[f'{tag}_RMS_all_kcal_A'] = round(s['rms_force_all_kcalmol_A'], 3)
            r[f'{tag}_RMS_solute_kcal_A'] = round(s['rms_force_solute_kcalmol_A'], 3)
        r['backbone_rmsd_A'] = round(mj['backbone_rmsd_A'], 3)
        r['mg_shell_final'] = ' | '.join(','.join(f'{n}:{dd}' for n, dd in sh) for sh in mj.get('mg_shell_final', []))
        if 'interface' in mj:
            i = mj['interface']
            r.update(if_contacts_start=i['heavy_atom_contacts_4.5A_start'], if_contacts_final=i['heavy_atom_contacts_4.5A_final'],
                     if_residues_start=i['interface_residues_start'], if_residues_final=i['interface_residues_final'],
                     if_CA_rmsd_A=round(i['interface_CA_rmsd_A'], 3))
    if len(r) > 1:
        rows.append(r)
keys = []
for r in rows:
    for k in r:
        if k not in keys:
            keys.append(k)
with open(out + '.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); [w.writerow(r) for r in rows]
json.dump(rows, open(out + '.json', 'w'), indent=1)
print(len(rows), 'systems ->', out + '.csv')
