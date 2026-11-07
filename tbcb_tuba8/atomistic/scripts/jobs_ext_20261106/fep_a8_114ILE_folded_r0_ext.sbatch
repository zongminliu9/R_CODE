#!/bin/bash
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 14:00:00
#SBATCH -J fep_a8_114ILE_folded_r0_ext
# S114I (TUBA8) folded r0 failed the leg half-split gate at 5000 iterations (0.74 > 0.27); extended to 10000
# iterations by the 2026-11-04 rule, as was done for I114S.  The leg was set up from the pre-2026-11-05 a8 build, so
# that build (superseded_20261105/build/a8) is passed again: a05 resumes from the checkpoint and the run directory's
# own wt_system.cif.gz / equil_state.npz; the build is only the side-chain donor library.
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
OLD=$SCRATCH/tbcb_atomistic/superseded_20261105/build/a8
R=$SCRATCH/tbcb_atomistic/runs/fep/a8_114ILE_folded/r0
[ -f $R/fep_ck_it5000.npz ] || cp -p $R/fep_ck.npz $R/fep_ck_it5000.npz
[ -f $R/fep_status_it5000.json ] || cp -p $R/fep_status.json $R/fep_status_it5000.json
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
date
python3 -W ignore scripts/a05_fep.py $OLD $SCRATCH/tbcb_atomistic/runs/fep/a8_114ILE_folded 0 --res 114 --to ILE --leg folded --ns 10.0 --budget 800
date
