#!/bin/bash
# 450 K unfolding analysis (a18) of the finished md450 replicas of one system into an interim directory.
# (2026-11-05: a18 failed in interim pass 5 because inputs/folding_units.json had never been copied to Sherlock.)
# Usage: sbatch -p normal -c 2 --mem=16G -t 0:40:00 -o <log> --wrap "bash -l scripts/run_a18_md450.sh <sid> <stamp>"
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
cd $HOME/tbcb_atomistic
SID=$1; O=$SCRATCH/tbcb_atomistic/interim/$2; mkdir -p $O
python3 -W ignore scripts/a18_unfold_analysis.py $SCRATCH/tbcb_atomistic/runs/md450/$SID $SCRATCH/tbcb_atomistic/build/$SID \
  $O/unfold_md450_$SID.json > $O/a18_md450_$SID.log 2>&1
echo "a18 exit $?" >> $O/a18_md450_$SID.log
