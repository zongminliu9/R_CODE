#!/bin/bash
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -J fep_e_a1a_cplx_451PHE_r2
# Rerun of TBCE-TUBA1A tail Y451F complex r2.  The first run (Sherlock 51130437, node sh03-12n12) lost its GPU
# ("NVML: GPU is lost", 2026-11-05 12:35) right after equilibration, hung without writing a checkpoint and timed out at
# 47 h.  The old run directory was moved to superseded_20261107/ and the leg starts again from equilibration; the time
# limit is 12 h (r0 / r1 took 3-5 h) so that a hang on a faulty node costs at most 12 h.  Same arguments otherwise.
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
B=$SCRATCH/tbcb_atomistic/build
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
date
python3 -W ignore scripts/a05_fep.py $B/e_a1a_cplx $SCRATCH/tbcb_atomistic/runs/fep/e_a1a_cplx_451PHE_folded 2 --res 451 --to PHE --leg folded --ns 5.0 --donor_lib $B/a1a/minimized.cif.gz --restrain "A:31,32,48,49,68,69|B:451" --budget 690
date
