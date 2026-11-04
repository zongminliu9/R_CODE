#!/usr/bin/env python3
"""Per-residue interaction-energy and force decomposition across the TBCB-alpha-tubulin
interface, for the TUBA8 and TUBA1A complexes.

What this computes, exactly. For a GIVEN structure the molecular-mechanics interaction between
two residues is an exact pairwise sum, with no sampling error:

    E_ij = sum over atom pairs [ 332.0637 q_a q_b / r_ab   +   4 eps_ab ( (s_ab/r)^12 - (s_ab/r)^6 ) ]

with Lorentz-Berthelot combination (s_ab = (s_a+s_b)/2, eps_ab = sqrt(eps_a eps_b)) and charges,
sigmas and epsilons read from the parametrised OpenMM system. Reported in kcal/mol.

This is an INTERACTION ENERGY DECOMPOSITION of one structure. It is not a free energy, it
contains no solvation screening (a separate GB term is reported for context), and it cannot be
converted to an affinity. Its use is to say WHICH residue pairs carry the interface and how that
differs between the two isoforms -- a question a single number for the whole interface cannot
answer.

Forces: the net force on each residue, and the residual force magnitude after minimisation, which
indicates local strain that the predicted structure carries.
"""
import argparse, glob, json, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import openmm as mm, openmm.app as app, openmm.unit as u
from pdbfixer import PDBFixer

COUL = 332.0637  # kcal/mol * A / e^2

ap = argparse.ArgumentParser()
ap.add_argument("--cif", required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--cutoff", type=float, default=12.0, help="A, for the pairwise sum")
ap.add_argument("--minimise", type=int, default=2000)
a = ap.parse_args()
os.makedirs("results/decomposition", exist_ok=True)

prep = f"data/models/prep/{a.tag}_decomp.pdb"
if not os.path.exists(prep):
    fx = PDBFixer(filename=a.cif)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.0)
    app.PDBFile.writeFile(fx.topology, fx.positions, open(prep,"w"), keepIds=True)
pdb = app.PDBFile(prep)
ff = app.ForceField("amber14-all.xml","implicit/gbn2.xml")
system = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                         nonbondedCutoff=1.8*u.nanometer, constraints=None)
atoms = list(pdb.topology.atoms())
N = len(atoms)
chains = sorted({at.residue.chain.id for at in atoms})
cnt = {c: sum(1 for at in atoms if at.residue.chain.id == c) for c in chains}
chT = min(cnt, key=cnt.get); chA = max(cnt, key=cnt.get)
print(f"[{a.tag}] {N} atoms; chain {chT} = TBCB ({cnt[chT]}), chain {chA} = alpha ({cnt[chA]})",
      flush=True)

integ = mm.VerletIntegrator(1*u.femtosecond)
sim = app.Simulation(pdb.topology, system, integ,
                     mm.Platform.getPlatformByName("CPU"), {"Threads":"1"})
sim.context.setPositions(pdb.positions)
e0 = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
f0 = sim.context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        u.kilocalorie_per_mole/u.angstrom)
print(f"[{a.tag}] as-predicted  E = {e0:.1f} kcal/mol, "
      f"max |F| = {np.linalg.norm(f0,axis=1).max():.1f} kcal/mol/A, "
      f"RMS |F| = {np.sqrt((f0**2).sum(1).mean()):.2f}", flush=True)
if a.minimise:
    sim.minimizeEnergy(maxIterations=a.minimise)
st = sim.context.getState(getPositions=True, getEnergy=True, getForces=True)
pos = st.getPositions(asNumpy=True).value_in_unit(u.angstrom)
frc = st.getForces(asNumpy=True).value_in_unit(u.kilocalorie_per_mole/u.angstrom)
e1 = st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
print(f"[{a.tag}] minimised     E = {e1:.1f} kcal/mol, "
      f"max |F| = {np.linalg.norm(frc,axis=1).max():.2f}, "
      f"RMS |F| = {np.sqrt((frc**2).sum(1).mean()):.3f}", flush=True)

nb = [system.getForce(i) for i in range(system.getNumForces())
      if isinstance(system.getForce(i), mm.NonbondedForce)][0]
q = np.zeros(N); sig = np.zeros(N); eps = np.zeros(N)
for i in range(N):
    c, s, e = nb.getParticleParameters(i)
    q[i] = c.value_in_unit(u.elementary_charge)
    sig[i] = s.value_in_unit(u.angstrom)
    eps[i] = e.value_in_unit(u.kilocalorie_per_mole)

idxT = np.array([i for i,at in enumerate(atoms) if at.residue.chain.id == chT])
idxA = np.array([i for i,at in enumerate(atoms) if at.residue.chain.id == chA])
resT = np.array([int(atoms[i].residue.id) for i in idxT])
resA = np.array([int(atoms[i].residue.id) for i in idxA])
nameT = {int(at.residue.id): at.residue.name for at in atoms if at.residue.chain.id == chT}
nameA = {int(at.residue.id): at.residue.name for at in atoms if at.residue.chain.id == chA}

pair = {}
CUT2 = a.cutoff**2
BLK = 400
for s in range(0, len(idxT), BLK):
    bi = idxT[s:s+BLK]
    d = pos[bi][:,None,:] - pos[idxA][None,:,:]
    r2 = (d*d).sum(-1)
    m = r2 < CUT2
    if not m.any(): continue
    ii, jj = np.where(m)
    r = np.sqrt(r2[ii,jj])
    ai = bi[ii]; aj = idxA[jj]
    ec = COUL*q[ai]*q[aj]/r
    sij = 0.5*(sig[ai]+sig[aj]); eij = np.sqrt(eps[ai]*eps[aj])
    x = np.divide(sij, r, out=np.zeros_like(r), where=r>0)
    x6 = x**6
    ev = 4.0*eij*(x6*x6 - x6)
    for k in range(len(r)):
        key = (int(resT[np.searchsorted(idxT, ai[k])]) if False else int(atoms[ai[k]].residue.id),
               int(atoms[aj[k]].residue.id))
        c, v = pair.get(key, (0.0, 0.0))
        pair[key] = (c + float(ec[k]), v + float(ev[k]))

rows = []
for (rt, ra), (c, v) in pair.items():
    rows.append(dict(tbcb_res=rt, tbcb_name=nameT.get(rt,"?"), alpha_res=ra,
                     alpha_name=nameA.get(ra,"?"), elec=round(c,4), vdw=round(v,4),
                     total=round(c+v,4)))
rows.sort(key=lambda r: r["total"])
tot_e = sum(r["elec"] for r in rows); tot_v = sum(r["vdw"] for r in rows)
print(f"[{a.tag}] pairwise interface sum over {len(rows)} residue pairs within {a.cutoff} A:")
print(f"[{a.tag}]   electrostatic {tot_e:10.2f}   van der Waals {tot_v:9.2f}   "
      f"total {tot_e+tot_v:10.2f} kcal/mol", flush=True)
print(f"[{a.tag}] strongest 12 residue pairs (kcal/mol):")
for r in rows[:12]:
    print(f"[{a.tag}]   TBCB {r['tbcb_name']}{r['tbcb_res']:<4d} -- "
          f"alpha {r['alpha_name']}{r['alpha_res']:<4d}  "
          f"elec {r['elec']:9.2f}  vdw {r['vdw']:7.2f}  total {r['total']:9.2f}")

# per-residue totals on the alpha side, and net force per residue
by_alpha = {}
for r in rows:
    k = r["alpha_res"]; c, v = by_alpha.get(k, (0.0,0.0))
    by_alpha[k] = (c + r["elec"], v + r["vdw"])
fres = {}
for i in idxA:
    k = int(atoms[i].residue.id); f = fres.get(k, np.zeros(3)); fres[k] = f + frc[i]
alpha_rows = [dict(alpha_res=k, alpha_name=nameA.get(k,"?"), elec=round(v[0],4),
                   vdw=round(v[1],4), total=round(v[0]+v[1],4),
                   net_force=round(float(np.linalg.norm(fres.get(k, np.zeros(3)))),4))
              for k, v in sorted(by_alpha.items())]
json.dump(dict(tag=a.tag, n_atoms=N, chain_tbcb=chT, chain_alpha=chA,
               E_as_predicted=e0, E_minimised=e1,
               maxF_as_predicted=float(np.linalg.norm(f0,axis=1).max()),
               rmsF_as_predicted=float(np.sqrt((f0**2).sum(1).mean())),
               maxF_minimised=float(np.linalg.norm(frc,axis=1).max()),
               rmsF_minimised=float(np.sqrt((frc**2).sum(1).mean())),
               cutoff_A=a.cutoff, n_pairs=len(rows),
               total_elec=tot_e, total_vdw=tot_v, total=tot_e+tot_v,
               pairs=rows, alpha_per_residue=alpha_rows),
          open(f"results/decomposition/{a.tag}.json","w"), indent=1)
print(f"[{a.tag}] wrote results/decomposition/{a.tag}.json", flush=True)
