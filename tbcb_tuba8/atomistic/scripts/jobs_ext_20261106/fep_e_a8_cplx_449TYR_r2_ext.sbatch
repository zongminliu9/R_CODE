#!/bin/bash
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -J fep_e_a8_cplx_449TYR_r2_ext
# Tail FEP replicate that failed the leg half-split gate at 5000 iterations (sixth interim pass, 2026-11-06); extended
# to 10000 iterations by the 2026-11-04 rule.  Same arguments as the original job; a05 resumes from the checkpoint.
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
B=$SCRATCH/tbcb_atomistic/build; F=$SCRATCH/tbcb_atomistic/runs/fep
R=$F/e_a8_cplx_449TYR_folded/r2
[ -f $R/fep_ck_it5000.npz ] || cp -p $R/fep_ck.npz $R/fep_ck_it5000.npz
[ -f $R/fep_status_it5000.json ] || cp -p $R/fep_status.json $R/fep_status_it5000.json
nvidia-smi --query-gpu=name,driver_version --format=csv,noheader
date
python3 -W ignore scripts/a05_fep.py $B/e_a8_cplx $F/e_a8_cplx_449TYR_folded 2 --res 449 --to TYR --leg folded --ns 10.0 --donor_lib $B/a1a/minimized.cif.gz --restrain "A:31,32,48,49,68,69|B:449" --budget 690
date
