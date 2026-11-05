"""Re-solvate a build's minimised solute in a larger box (for the unbinding PMF), then minimise with force record.
Usage: a16_resolvate.py <src_build> <new_build> --pad 2.2"""
import os, sys, json, gzip, argparse, subprocess
import numpy as np
import openmm as mm
from openmm import app, unit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
ap = argparse.ArgumentParser(); ap.add_argument('src'); ap.add_argument('dst'); ap.add_argument('--pad', type=float, default=2.2)
a = ap.parse_args()
top, pos = C.load_cif(os.path.join(a.src, 'minimized.cif.gz'))
m = app.Modeller(top, pos)
m.delete([r for r in m.topology.residues() if r.name in C.SOLVENT])
ff = C.forcefield()
m.addSolvent(ff, model='tip3p', padding=a.pad * C.NM, boxShape='dodecahedron', ionicStrength=0.15 * unit.molar, neutralize=True)
os.makedirs(a.dst, exist_ok=True)
C.write_cif(os.path.join(a.dst, 'solvated.cif.gz'), m.topology, m.positions)
json.dump(dict(src=a.src, pad_nm=a.pad, n_atoms=m.topology.getNumAtoms()), open(os.path.join(a.dst, 'build.json'), 'w'), indent=1)
subprocess.run([sys.executable, '-W', 'ignore', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a02_minimize.py'), a.dst], check=True)
print(a.dst, m.topology.getNumAtoms())
