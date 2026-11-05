"""HADDOCK3 docking of TBCB (UBL + linker, residues 1-133) onto an alpha-tubulin body (1-437 with GTP.Mg2+), for the
ensemble-docking step (native structure, then representative intermediates from REST2/MSM).

Restrained mode: ambiguous restraints between the alpha residues that both AF3 predictions place in the interface
(robust in all seeds, alpha <= 437; R_CODE results/af3/interface_aligned.json) and the TBCB residues within 5 A of
them in the AF3 model.  The divergent (TUBA8-specific) positions are never used as restraints.
Ab initio mode: no interface information (centre-of-mass restraints + random surface AIRs), as the control that
tells whether the restraints drive the answer.
Scored against the AF3 pose (caprieval: i-RMSD, l-RMSD, fnat) so we can say whether docking finds the AF3 site.
HADDOCK scores are reported as scores only, never as energies.

Usage: a14_dock.py prep <alpha_source> <out_dir> [--frame cif.gz/pdb]   (alpha_source: a build dir or a pdb)
       a14_dock.py cfg  <out_dir> --mode restrained|abinitio [--sampling 1000] [--ncores 16]
"""
import os, sys, json, argparse, subprocess
import numpy as np
import gemmi

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
HD = '/scratch/m000251/zongminl/mm/envs/hd3/bin'
ap = argparse.ArgumentParser()
ap.add_argument('mode'); ap.add_argument('a1'); ap.add_argument('a2', nargs='?')
ap.add_argument('--dock', default='restrained'); ap.add_argument('--sampling', type=int, default=1000)
ap.add_argument('--ncores', type=int, default=16); ap.add_argument('--af3', default='fold_tbcb_tuba1a')
a = ap.parse_args()
SHARED = None


def shared_alpha_interface():
    j = json.load(open(os.path.join(ROOT, 'inputs', 'interface_aligned.json')))
    s = set(j['TUBA8']['robust_canon']) & set(j['TUBA1A']['robust_canon'])
    return sorted(p for p in s if p <= 437)


def write_alpha(src, out):
    """protein 1-437 heavy atoms + GTP heavy atoms + Mg as MG2 (charged Mg2+ in HADDOCK's ion.top), chain A."""
    st = gemmi.read_structure(src); st.setup_entities(); st.remove_hydrogens(); st.remove_waters()
    o = gemmi.Structure(); o.add_model(gemmi.Model('1')); ch = gemmi.Chain('A')
    for c in st[0]:
        for r in c:
            if r.name in ('NA', 'CL', 'HOH', 'MG'):
                # Mg2+ is dropped for docking: CNS torsion-angle dynamics (flexref) aborts on an isolated atom
                # ('%TORSION:TOPOLOGY error ... isolated (not covalently bonded) atoms').  GTP stays.
                continue
            tab = gemmi.find_tabulated_residue(r.name)
            if r.name in ('GTP', 'MG') or (tab is not None and tab.is_amino_acid()):
                nr = gemmi.Residue(); nr.name = 'MG2' if r.name == 'MG' else r.name; nr.seqid = r.seqid
                if r.name == 'GTP': nr.seqid = gemmi.SeqId(501, ' ')
                if r.name == 'MG': nr.seqid = gemmi.SeqId(502, ' ')
                nr.het_flag = 'H' if r.name in ('GTP', 'MG') else 'A'
                for at in r:
                    if at.name == 'OXT':
                        continue
                    nr.add_atom(at)
                ch.add_residue(nr)
    o[0].add_chain(ch); o.setup_entities()
    o.write_pdb(out)


if a.mode == 'prep':
    out = a.a2; os.makedirs(out, exist_ok=True)
    src = a.a1 if a.a1.endswith(('.pdb', '.cif', '.gz')) else os.path.join(a.a1, 'start_noH.pdb')
    write_alpha(src, os.path.join(out, 'alpha.pdb'))
    af = gemmi.read_structure(os.path.join(ROOT, 'inputs', 'af3', a.af3, f'{a.af3}_model_0.cif')); af.remove_hydrogens()
    tb = gemmi.Structure(); tb.add_model(gemmi.Model('1')); cb = gemmi.Chain('B')
    for r in af[0]['A'].get_polymer():
        if r.seqid.num <= 133:
            cb.add_residue(r)
    tb[0].add_chain(cb); tb.setup_entities(); tb.write_pdb(os.path.join(out, 'tbcb.pdb'))
    # reference complex (AF3 alpha body + TBCB 1-133) for caprieval
    ref = gemmi.Structure(); ref.add_model(gemmi.Model('1')); ca = gemmi.Chain('A')
    for r in af[0]['B'].get_polymer():
        if r.seqid.num <= 437:
            ca.add_residue(r)
    ref[0].add_chain(ca); ref[0].add_chain(cb); ref.setup_entities(); ref.write_pdb(os.path.join(out, 'ref.pdb'))
    act_a = shared_alpha_interface()
    ns = gemmi.NeighborSearch(af[0], af.cell, 6).populate()
    act_b = set()
    for r in af[0]['B'].get_polymer():
        if r.seqid.num in act_a:
            for at in r:
                for mk in ns.find_atoms(at.pos, '\0', radius=5.0):
                    cra = mk.to_cra(af[0])
                    if cra.chain.name == 'A' and cra.residue.seqid.num <= 133:
                        act_b.add(cra.residue.seqid.num)
    act_b = sorted(act_b)
    open(os.path.join(out, 'alpha.act'), 'w').write(' '.join(map(str, act_a)) + '\n\n')
    open(os.path.join(out, 'tbcb.act'), 'w').write(' '.join(map(str, act_b)) + '\n\n')
    with open(os.path.join(out, 'ambig.tbl'), 'w') as fh:
        subprocess.run([f'{HD}/haddock3-restraints', 'active_passive_to_ambig', os.path.join(out, 'alpha.act'),
                        os.path.join(out, 'tbcb.act'), '--segid-one', 'A', '--segid-two', 'B'], stdout=fh, check=True)
    json.dump(dict(alpha_active=act_a, tbcb_active=act_b, src=src), open(os.path.join(out, 'prep.json'), 'w'), indent=1)
    print(len(act_a), 'alpha active;', len(act_b), 'TBCB active')
elif a.mode == 'cfg':
    out = a.a1
    base = f"""run_dir = "run_{a.dock}"
ncores = {a.ncores}
mode = "local"
molecules = ["alpha.pdb", "tbcb.pdb"]

[topoaa]
autohis = true

"""
    if a.dock == 'restrained':
        body = f"""[rigidbody]
ambig_fname = "ambig.tbl"
sampling = {a.sampling}

[caprieval]
reference_fname = "ref.pdb"

[seletop]
select = 200

[flexref]
ambig_fname = "ambig.tbl"

[emref]
ambig_fname = "ambig.tbl"

[clustfcc]

[seletopclusts]
top_models = 4

[caprieval]
reference_fname = "ref.pdb"
"""
    else:
        body = f"""[rigidbody]
ranair = true
cmrest = true
sampling = {2 * a.sampling}

[caprieval]
reference_fname = "ref.pdb"

[seletop]
select = 200

[flexref]
contactairs = true

[emref]
contactairs = true

[clustfcc]

[seletopclusts]
top_models = 4

[caprieval]
reference_fname = "ref.pdb"
"""
    open(os.path.join(out, f'{a.dock}.cfg'), 'w').write(base + body)
    print('wrote', os.path.join(out, f'{a.dock}.cfg'))


def analyse(out_dir, mode):
    """Cluster-level results of the last caprieval step: HADDOCK score (a score, not an energy), i-RMSD/l-RMSD/fnat/
    DockQ against the AF3 pose, cluster size.  Writes dock_<mode>.json."""
    import csv, glob
    rd = os.path.join(out_dir, f'run_{mode}')
    cap = sorted(glob.glob(os.path.join(rd, '*_caprieval')))
    if not cap:
        print('no caprieval in', rd); return None
    last = cap[-1]
    rows = []
    with open(os.path.join(last, 'capri_clt.tsv')) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'):
            rows.append(r)
    ss = []
    with open(os.path.join(last, 'capri_ss.tsv')) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith('#')), delimiter='\t'):
            ss.append(r)
    res = dict(mode=mode, n_clusters=len(rows), clusters=rows[:10], best_models=ss[:10])
    json.dump(res, open(os.path.join(out_dir, f'dock_{mode}.json'), 'w'), indent=1)
    for r in rows[:5]:
        print(mode, {k: r.get(k) for k in ('cluster_rank', 'n', 'score', 'irmsd', 'fnat', 'lrmsd', 'dockq')})
    return res


if a.mode == 'analyse':
    analyse(a.a1, a.dock)
