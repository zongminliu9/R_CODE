"""MM/GBSA and MM/PBSA (AmberTools MMPBSA.py, single-trajectory protocol) on the explicit-solvent MD of a complex,
with per-residue decomposition.  Topologies are converted from the OpenMM system with ParmEd (ff14SB + Meagher GTP
+ IOD Mg2+ + TIP3P parameters exactly as simulated), solute only; mbondi2 radii.  The four Mg-shell waters are kept
with the tubulin they belong to.

These numbers rank poses and isoforms within one binding mode; they are NOT affinities and are never converted to Kd.
The uncertainty quoted is the spread across independent MD replicas, not the per-frame standard error.

Usage: a12_mmpbsa.py prep <build_dir> <out_dir> --receptor A --ligand B
       a12_mmpbsa.py run  <out_dir> <md_run_dir/rK> [--start_ns 10] [--n_gb 250] [--n_pb 50] [--decomp]
"""
import os, sys, json, glob, argparse, subprocess, shutil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('mode'); ap.add_argument('a1'); ap.add_argument('a2')
ap.add_argument('--receptor', default='A'); ap.add_argument('--ligand', default='B')
ap.add_argument('--start_ns', type=float, default=10.0); ap.add_argument('--n_gb', type=int, default=250)
ap.add_argument('--n_pb', type=int, default=50); ap.add_argument('--decomp', action='store_true')
ap.add_argument('--dt_ps', type=float, default=20.0)
a = ap.parse_args()
AMB = '/scratch/m000251/zongminl/mm/envs/amb'
if a.mode == 'prep':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import common as C
    import parmed
    import openmm as mm
    from openmm import app, unit
    os.makedirs(a.a2, exist_ok=True)
    top, pos = C.load_cif(os.path.join(a.a1, 'minimized.cif.gz'))
    shell = {b for _, b in C.mg_shell_pairs(top, pos)}
    shell_res = {list(top.atoms())[b].residue.index for b in shell}
    keep = [x.index for x in top.atoms() if x.residue.name not in C.SOLVENT or x.residue.index in shell_res]
    m = app.Modeller(top, pos)
    ks = set(keep)
    m.delete([x for x in m.topology.atoms() if x.index not in ks])
    ff = C.forcefield()
    s = ff.createSystem(m.topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    st = parmed.openmm.load_topology(m.topology, s, xyz=m.positions)
    action = parmed.tools.changeRadii(st, 'mbondi2'); action.execute()
    # chain membership: receptor chains, ligand chains (+ GTP/Mg/shell waters, which follow the tubulin chain)
    rec = set(a.receptor.split(',')); lig = set(a.ligand.split(','))
    def chain_of(r):
        return r.chain
    nres = len(st.residues)
    ridx_rec = [i for i, r in enumerate(st.residues) if r.chain in rec]
    ridx_lig = [i for i in range(nres) if i not in set(ridx_rec)]
    st.save(os.path.join(a.a2, 'complex.prmtop'), overwrite=True)
    st.save(os.path.join(a.a2, 'complex.pdb'), overwrite=True)
    def sub(idx, fn):
        mask = '@' + ','.join(str(at.idx + 1) for i in idx for at in st.residues[i].atoms)
        sst = st[mask]; sst.save(os.path.join(a.a2, fn), overwrite=True)
        return len(sst.atoms)
    nr = sub(ridx_rec, 'receptor.prmtop'); nl = sub(ridx_lig, 'ligand.prmtop')
    # atom order of the MD solute selection must equal this complex (solute.pdb written by a03_md.py)
    info = dict(n_complex=len(st.atoms), n_receptor=nr, n_ligand=nl, receptor_chains=sorted(rec),
                receptor_residue_range=[ridx_rec[0] + 1, ridx_rec[-1] + 1] if ridx_rec else None,
                ligand_residue_range=[ridx_lig[0] + 1, ridx_lig[-1] + 1] if ridx_lig else None)
    # residue masks for MMPBSA (1-based residue numbers in the stripped complex)
    def ranges(ix):
        ix = [i + 1 for i in ix]; out = []; s0 = ix[0]; p = ix[0]
        for v in ix[1:]:
            if v != p + 1:
                out.append(f'{s0}-{p}'); s0 = v
            p = v
        out.append(f'{s0}-{p}'); return ','.join(out)
    info['receptor_mask'] = ':' + ranges(ridx_rec); info['ligand_mask'] = ':' + ranges(ridx_lig)
    json.dump(info, open(os.path.join(a.a2, 'prep.json'), 'w'), indent=1)
    print(json.dumps(info))
else:
    out, rdir = a.a1, a.a2
    info = json.load(open(os.path.join(out, 'prep.json')))
    segs = sorted(glob.glob(os.path.join(rdir, 'seg_*.dcd')))
    tag = os.path.basename(rdir.rstrip('/'))
    W = os.path.join(out, tag); os.makedirs(W, exist_ok=True)
    # concatenate segments with cpptraj, keeping only frames after start_ns
    first = int(a.start_ns * 1000 / a.dt_ps) + 1
    cin = [f'parm {os.path.join(out, "complex.prmtop")}']
    for sgm in segs:
        cin.append(f'trajin {sgm}')
    cin += [f'trajout {W}/prod.nc netcdf start {first}', 'run']
    open(os.path.join(W, 'cat.in'), 'w').write('\n'.join(cin) + '\n')
    subprocess.run([f'{AMB}/bin/cpptraj', '-i', os.path.join(W, 'cat.in')], check=True, stdout=open(os.path.join(W, 'cat.log'), 'w'))
    nfr = int(subprocess.run([f'{AMB}/bin/cpptraj', '-p', os.path.join(out, 'complex.prmtop'), '-y', f'{W}/prod.nc', '-tl'],
                             capture_output=True, text=True).stdout.split()[-1])
    for kind, nwant in (('gb', a.n_gb), ('pb', a.n_pb)):
        itv = max(1, nfr // max(1, nwant))
        block = ('&gb\n igb=5, saltcon=0.150,\n/\n' if kind == 'gb' else '&pb\n istrng=0.150, inp=2, radiopt=0, fillratio=4.0,\n/\n')
        dec = '&decomp\n idecomp=2, dec_verbose=0,\n/\n' if (a.decomp and kind == 'gb') else ''
        txt = (f"&general\n startframe=1, endframe={nfr}, interval={itv}, keep_files=0, verbose=1,\n"
               f" receptor_mask='{info['receptor_mask']}', ligand_mask='{info['ligand_mask']}',\n/\n{block}{dec}")
        inp = os.path.join(W, f'{kind}.in'); open(inp, 'w').write(txt)
        cmd = [f'{AMB}/bin/MMPBSA.py', '-O', '-i', inp, '-o', os.path.join(W, f'{kind}_FINAL.dat'),
               '-eo', os.path.join(W, f'{kind}_frames.csv'), '-cp', os.path.join(out, 'complex.prmtop'),
               '-rp', os.path.join(out, 'receptor.prmtop'), '-lp', os.path.join(out, 'ligand.prmtop'), '-y', f'{W}/prod.nc']
        if dec:
            cmd += ['-do', os.path.join(W, 'decomp.dat'), '-deo', os.path.join(W, 'decomp_frames.csv')]
        subprocess.run(cmd, cwd=W, check=True, stdout=open(os.path.join(W, f'{kind}.log'), 'w'), stderr=subprocess.STDOUT)
    print('done', W, nfr, 'frames')
