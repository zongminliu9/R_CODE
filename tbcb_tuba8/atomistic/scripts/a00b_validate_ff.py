"""Check that ISOLDE's ff14SB conversion (used here because the GTP template shares its unprefixed Amber types)
gives the same protein energy, term by term, as OpenMM's own amber14/protein.ff14SB.xml.  Protein only, vacuum,
no cutoff, at the minimised coordinates of one build."""
import os, sys, gzip, json
import numpy as np
import openmm as mm
from openmm import app, unit
FF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ff')
app.Topology.loadBondDefinitions(os.path.join(FF, 'gtp_bonds.xml'))
with gzip.open(sys.argv[1], 'rt') as fh:
    px = app.PDBxFile(fh)
m = app.Modeller(px.topology, px.positions)
m.delete([r for r in m.topology.residues() if r.name in ('HOH', 'NA', 'CL', 'MG', 'GTP')])
out = {}
for tag, files in (('isolde', [os.path.join(FF, 'amberff14SB.xml')]), ('openmm', [os.path.join(FF, 'openmm_protein.ff14SB.xml')])):
    ff = app.ForceField(*files)
    s = ff.createSystem(m.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    for i, f in enumerate(s.getForces()): f.setForceGroup(i)
    c = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName('Reference'))
    c.setPositions(m.positions)
    out[tag] = {type(f).__name__: c.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
                for i, f in enumerate(s.getForces()) if type(f).__name__ != 'CMMotionRemover'}
    out[tag]['total'] = c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    q = [f for f in s.getForces() if isinstance(f, mm.NonbondedForce)][0]
    out[tag]['net_charge'] = sum(q.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(s.getNumParticles()))
for k in out['isolde']:
    print(f"{k:24s} isolde {out['isolde'][k]:14.4f}  openmm {out['openmm'].get(k, float('nan')):14.4f}  diff {out['isolde'][k]-out['openmm'].get(k, float('nan')):.2e}")
json.dump(out, open(os.path.join(os.path.dirname(sys.argv[1]), 'ff_validation.json'), 'w'), indent=1)
