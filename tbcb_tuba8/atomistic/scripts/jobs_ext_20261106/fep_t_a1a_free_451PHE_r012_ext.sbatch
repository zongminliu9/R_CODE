#!/bin/bash
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 14:00:00
#SBATCH -J fep_t_a1a_free_451PHE_r012_ext
# Free TUBA1A tail Y451F legs r0-r2 failed the leg half-split gate at 5000 iterations (0.62 / 0.21 / 0.21 > 0.114);
# each is extended to 10000 iterations by the 2026-11-04 rule, one after the other in one GPU job.
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
B=$SCRATCH/tbcb_atomistic/build; F=$SCRATCH/tbcb_atomistic/runs/fep/t_a1a_free_451PHE_folded
for r in 0 1 2; do
  [ -f $F/r$r/fep_ck_it5000.npz ] || cp -p $F/r$r/fep_ck.npz $F/r$r/fep_ck_it5000.npz
  [ -f $F/r$r/fep_status_it5000.json ] || cp -p $F/r$r/fep_status.json $F/r$r/fep_status_it5000.json
done
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
date
T0=$(date +%s)
for r in 0 1 2; do
  LEFT=$(( 800 - ($(date +%s) - T0) / 60 ))
  [ $LEFT -gt 30 ] && python3 -W ignore scripts/a05_fep.py $B/t_a1a_free $F $r --res 451 --to PHE --leg folded --ns 10.0 --donor_lib $B/a1a/minimized.cif.gz --budget $LEFT
done
date
