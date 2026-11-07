#!/bin/bash
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -J fep_a8_114ILE_unfolded_r012_ext
# S114I (TUBA8) unfolded r0, r1, r2 failed the leg half-split gate at 5000 iterations (0.24 / 0.33 / 0.51 > 0.239);
# each is extended to 10000 iterations by the 2026-11-04 rule, one after the other in one GPU job.  r0 and r1 were set
# up from the pre-2026-11-05 a8 build, r2 from the rebuilt one; each is resumed with the build it was set up from
# (only the side-chain donor library is read from it).
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
OLD=$SCRATCH/tbcb_atomistic/superseded_20261105/build/a8
NEW=$SCRATCH/tbcb_atomistic/build/a8
U=$SCRATCH/tbcb_atomistic/runs/fep/a8_114ILE_unfolded
for r in 0 1 2; do
  [ -f $U/r$r/fep_ck_it5000.npz ] || cp -p $U/r$r/fep_ck.npz $U/r$r/fep_ck_it5000.npz
  [ -f $U/r$r/fep_status_it5000.json ] || cp -p $U/r$r/fep_status.json $U/r$r/fep_status_it5000.json
done
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
date
T0=$(date +%s)
for r in 0 1 2; do
  B=$OLD; [ $r -eq 2 ] && B=$NEW
  LEFT=$(( 700 - ($(date +%s) - T0) / 60 ))
  [ $LEFT -gt 30 ] && python3 -W ignore scripts/a05_fep.py $B $U $r --res 114 --to ILE --leg unfolded --ns 10.0 --budget $LEFT
done
date
