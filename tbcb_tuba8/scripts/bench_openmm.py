#!/usr/bin/env python3
"""Measure real OpenMM CPU throughput for system sizes relevant to this project."""
import time, numpy as np, openmm as mm, openmm.app as app, openmm.unit as u

def build_ala(n, implicit):
    """n-residue poly-ALA as a size surrogate; returns (system, context-ready sim)."""
    top = app.Topology(); ch = top.addChain(); positions=[]
    # crude extended backbone+CB
    for i in range(n):
        r = top.addResidue("ALA", ch)
        names=[("N","N"),("CA","C"),("CB","C"),("C","C"),("O","O")]
        atoms=[]
        for nm,el in names:
            atoms.append(top.addAtom(nm, app.Element.getBySymbol(el), r))
        top.addBond(atoms[0],atoms[1]); top.addBond(atoms[1],atoms[2])
        top.addBond(atoms[1],atoms[3]); top.addBond(atoms[3],atoms[4])
        if i>0: top.addBond(prevC, atoms[0])
        prevC = atoms[3]
        x=i*0.36
        positions += [[x,0,0],[x+0.12,0.05,0],[x+0.12,0.19,0],[x+0.25,0.0,0],[x+0.25,-0.12,0]]
    return top, np.array(positions)

print(f"{'system':34s} {'atoms':>7s} {'ns/day':>9s} {'days/100ns':>11s}")
ff_imp = app.ForceField('amber14-all.xml','implicit/gbn2.xml')
for n,label in [(450,'alpha-tubulin monomer, GBn2'),(90,'CAP-Gly domain, GBn2'),(100,'CAP-Gly+tail peptide, GBn2')]:
    top,pos = build_ala(n, True)
    try:
        sysobj = ff_imp.createSystem(top, nonbondedMethod=app.CutoffNonPeriodic,
                    nonbondedCutoff=1.8*u.nanometer, constraints=app.HBonds, hydrogenMass=None,
                    soluteDielectric=1.0, solventDielectric=78.5)
    except Exception as e:
        print(f"  {label}: build fail {e}"); continue
    integ = mm.LangevinMiddleIntegrator(310*u.kelvin, 1/u.picosecond, 0.004*u.picoseconds)
    plat = mm.Platform.getPlatformByName('CPU')
    sim = app.Simulation(top, sysobj, integ, plat, {'Threads':'4'})
    sim.context.setPositions(pos)
    sim.minimizeEnergy(maxIterations=200)
    sim.step(100)                      # warm-up
    t0=time.time(); nsteps=1500; sim.step(nsteps); dt=time.time()-t0
    ns = nsteps*0.004/1000.0; nsday = ns/dt*86400
    print(f"{label:34s} {sysobj.getNumParticles():7d} {nsday:9.1f} {100/nsday:11.1f}")
