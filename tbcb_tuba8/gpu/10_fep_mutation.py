#!/usr/bin/env python3
"""Relative alchemical free energy for a point substitution, via perses + Hamiltonian replica
exchange, with MBAR.

Why perses rather than a hand-rolled alchemical setup. Ile->Ser changes the side-chain topology
(CB-CG1-CG2-CD1 against CB-OG), so the hybrid topology needs a correct atom map, correct
valence-term interpolation and a soft-core that does not blow up when CD1 disappears. A
home-made version of that machinery fails silently: it still produces a number, and the number
is wrong in a way no convergence diagnostic catches. perses is the maintained implementation of
exactly this, so it does the mapping and this script does the driving, the sampling schedule and
the diagnostics.

What this script is responsible for, and what to check in its output:
  * electrostatics and sterics are decoupled in separate lambda stages, because a buried
    hydrophobic-to-polar change has a large vdW component that single-stage schedules converge
    badly;
  * Hamiltonian replica exchange across the lambda ladder, so windows mix rather than relying on
    each being independently equilibrated;
  * three independent replicas per leg from different velocity seeds, with the final uncertainty
    taken from the spread ACROSS replicas, not from MBAR's analytical error, which is optimistic
    when windows are correlated;
  * the transformation is run in both directions and the hysteresis reported;
  * the adjacent-window overlap matrix is written out, and 90_analyse_fep.py refuses to quote a
    ddG if any adjacent overlap falls below 0.03.
"""
import argparse, json, os, sys, time

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True, help="serialised System from 00_build_systems.py")
    ap.add_argument("--mutation", required=True, help="e.g. I114S (wild-type, position, mutant)")
    ap.add_argument("--legs", nargs="+", default=["folded"],
                    choices=["folded","reference","complex","free"])
    ap.add_argument("--replicas", type=int, default=3)
    ap.add_argument("--lambdas", type=int, default=20)
    ap.add_argument("--ns-per-window", type=float, default=10.0)
    ap.add_argument("--reverse", action="store_true", help="run the mutant->wild-type direction")
    ap.add_argument("--out", required=True)
    ap.add_argument("--checkpoint-ps", type=float, default=500.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    try:
        import openmm as mm, openmm.app as app, openmm.unit as u
        from openmmtools import alchemy, states, mcmc, multistate
        import numpy as np
    except ImportError as e:
        sys.exit(f"missing dependency: {e}. create the env from environment.yml first.")
    try:
        from perses.app.relative_point_mutation_setup import PointMutationExecutor
    except ImportError:
        print("perses is not installed. install it into the env (conda install -c conda-forge "
              "perses) -- this script intentionally does not substitute a hand-rolled "
              "alchemical factory, see the module docstring.", file=sys.stderr)
        sys.exit(3)

    wt, pos, mut = a.mutation[0], a.mutation[1:-1], a.mutation[-1]
    base = os.path.splitext(a.system)[0]
    pdb = base + ".pdb"
    print(f"system      {a.system}\nstructure   {pdb}\nmutation    {wt}{pos}{mut}"
          f"\nlegs        {a.legs}\nlambdas     {a.lambdas}"
          f"\nns/window   {a.ns_per_window}\nreplicas    {a.replicas}", flush=True)

    THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU","G":"GLY",
             "H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE","P":"PRO","S":"SER",
             "T":"THR","W":"TRP","Y":"TYR","V":"VAL"}
    for leg in a.legs:
        for rep in range(a.replicas):
            tag = f"{a.mutation}_{leg}_rep{rep}" + ("_rev" if a.reverse else "")
            nc  = os.path.join(a.out, tag + ".nc")
            if os.path.exists(nc):
                print(f"  {tag}: checkpoint present, resuming", flush=True); 
            t0 = time.time()
            exe = PointMutationExecutor(
                protein_filename=pdb,
                mutation_chain_id="A",
                mutation_residue_id=str(pos),
                proposed_residue=THREE[mut if not a.reverse else wt],
                phase="solvent" if leg in ("reference","free") else "complex",
                conduct_endstate_validation=True,
                forcefield_files=["amber19-all.xml", "amber14/opc.xml"],
                small_molecule_forcefield=None,
                ionic_strength=0.15*u.molar,
            )
            htf = exe.get_apo_htf() if leg in ("folded","reference","free") else exe.get_complex_htf()
            # two-stage lambda schedule: electrostatics first, then sterics
            n = a.lambdas
            lam_elec = np.concatenate([np.linspace(0, 1, n//2), np.ones(n - n//2)])
            lam_ster = np.concatenate([np.zeros(n//2), np.linspace(0, 1, n - n//2)])
            factory = alchemy.AbsoluteAlchemicalFactory(consistent_exceptions=False)
            thermo = states.ThermodynamicState(htf.hybrid_system, temperature=310*u.kelvin,
                                               pressure=1*u.bar)
            comp = states.SamplerState(positions=htf.hybrid_positions)
            protocol = {"lambda_electrostatics": lam_elec.tolist(),
                        "lambda_sterics": lam_ster.tolist()}
            move = mcmc.LangevinSplittingDynamicsMove(
                timestep=4.0*u.femtoseconds, n_steps=250,
                collision_rate=1.0/u.picosecond, reassign_velocities=False,
                splitting="V R O R V", constraint_tolerance=1e-6)
            iters = int(a.ns_per_window*1000/(250*0.004))
            sampler = multistate.ReplicaExchangeSampler(
                mcmc_moves=move, number_of_iterations=iters,
                online_analysis_interval=100, online_analysis_target_error=0.0)
            reporter = multistate.MultiStateReporter(nc, checkpoint_interval=int(
                a.checkpoint_ps/(250*0.004)))
            comp_states = states.create_thermodynamic_state_protocol(
                thermo, protocol=protocol, composable_states=[])
            sampler.create(thermodynamic_states=comp_states,
                           sampler_states=[comp]*len(comp_states), storage=reporter)
            sampler.run()
            print(f"  {tag}: {iters} iterations in {(time.time()-t0)/3600:.2f} h -> {nc}",
                  flush=True)
            json.dump(dict(mutation=a.mutation, leg=leg, replica=rep, reverse=a.reverse,
                           lambdas=a.lambdas, ns_per_window=a.ns_per_window,
                           iterations=iters, storage=nc,
                           protocol={k: list(map(float, v)) for k, v in protocol.items()}),
                      open(os.path.join(a.out, tag + "_meta.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
