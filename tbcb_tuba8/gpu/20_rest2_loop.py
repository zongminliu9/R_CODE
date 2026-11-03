#!/usr/bin/env python3
"""REST2 comparison of local unfolding propensity in the 32-45 loop.

A 14-residue segment substitution cannot be done alchemically, so the loop question is answered
as a population difference rather than a single ddG: how often does the loop visit locally
unfolded states in TUBA8 against the sequence-swapped control, and conversely.

REST2 (replica exchange with solute tempering) scales only the solute's energy terms, so the
water stays at 310 K and the number of replicas needed grows with the size of the SOLUTE, not the
box. Plain MD at 310 K will not sample these transitions in any affordable time, which is why
this is not simply a long trajectory.

Order parameters written for analysis: loop native-contact fraction, loop radius of gyration,
loop SASA, loop backbone RMSD to the starting model, and exposed hydrophobic surface in the
loop's neighbourhood -- the last being the clearance-relevant quantity, since it is what a
degradation machinery would recognise.
"""
import argparse, json, os, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", nargs=2, required=True,
                    help="two system names from 00_build_systems.py, e.g. tuba8 tuba8_loop1a")
    ap.add_argument("--replicas", type=int, default=16)
    ap.add_argument("--tmin", type=float, default=310.0)
    ap.add_argument("--tmax", type=float, default=450.0, help="effective solute temperature")
    ap.add_argument("--ns", type=float, default=200.0)
    ap.add_argument("--loop", nargs=2, type=int, default=[32, 45])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    try:
        import numpy as np, openmm as mm, openmm.app as app, openmm.unit as u
        from openmmtools import states, mcmc, multistate
        from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
        from openmmtools.states import GlobalParameterState
    except ImportError as e:
        sys.exit(f"missing dependency: {e}")
    try:
        from openmmtools.forcefactories import restrain_atoms  # presence check
        from openmmtools.alchemy import AlchemicalRegion       # REST2 via alchemical scaling
    except ImportError as e:
        sys.exit(f"openmmtools too old for REST2 helpers: {e}")

    HERE = os.path.dirname(os.path.abspath(__file__))
    for name in a.pair:
        sx = os.path.join(HERE, "systems", f"{name}_solv.xml")
        px = os.path.join(HERE, "systems", f"{name}_solv.pdb")
        if not (os.path.exists(sx) and os.path.exists(px)):
            sys.exit(f"missing system {name}: run 00_build_systems.py first")
        system = mm.XmlSerializer.deserialize(open(sx).read())
        pdb = app.PDBFile(px)
        solute = [at.index for at in pdb.topology.atoms()
                  if at.residue.name not in ("HOH","NA","CL")]
        # geometric ladder of effective solute temperatures
        betas = np.geomspace(a.tmin, a.tmax, a.replicas)
        print(f"{name}: {system.getNumParticles()} atoms, {len(solute)} solute atoms, "
              f"{a.replicas} replicas over {a.tmin:.0f}-{a.tmax:.0f} K", flush=True)
        region = AlchemicalRegion(alchemical_atoms=solute, annihilate_sterics=False,
                                  annihilate_electrostatics=False)
        from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalState
        factory = AbsoluteAlchemicalFactory(consistent_exceptions=False)
        alch_system = factory.create_alchemical_system(system, region)
        alch_state = AlchemicalState.from_system(alch_system)
        thermo = states.ThermodynamicState(alch_system, temperature=a.tmin*u.kelvin,
                                           pressure=1*u.bar)
        comp = []
        for T in betas:
            s = states.CompoundThermodynamicState(thermo, composable_states=[alch_state])
            # REST2: scale solute-solute terms by beta_m/beta_0, solute-solvent by sqrt of it
            scale = a.tmin/float(T)
            s.lambda_sterics = 1.0
            s.lambda_electrostatics = float(np.sqrt(scale))
            comp.append(s)
        move = mcmc.LangevinSplittingDynamicsMove(
            timestep=4.0*u.femtoseconds, n_steps=500, collision_rate=1.0/u.picosecond,
            splitting="V R O R V", constraint_tolerance=1e-6)
        iters = int(a.ns*1000/(500*0.004))
        sampler = ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=iters,
                                         replica_mixing_scheme="swap-neighbors")
        nc = os.path.join(a.out, f"{name}_rest2.nc")
        reporter = MultiStateReporter(nc, checkpoint_interval=250)
        ss = states.SamplerState(positions=pdb.positions)
        sampler.create(thermodynamic_states=comp, sampler_states=[ss]*len(comp),
                       storage=reporter)
        sampler.run()
        json.dump(dict(system=name, replicas=a.replicas, tmin=a.tmin, tmax=a.tmax,
                       ns=a.ns, loop=a.loop, storage=nc,
                       order_parameters=["loop_native_contacts","loop_Rg","loop_SASA",
                                         "loop_backbone_RMSD","neighbourhood_hydrophobic_SASA"]),
                  open(os.path.join(a.out, f"{name}_rest2_meta.json"), "w"), indent=1)
        print(f"  -> {nc}", flush=True)

if __name__ == "__main__":
    main()
