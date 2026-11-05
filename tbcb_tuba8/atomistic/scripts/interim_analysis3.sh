#!/bin/bash
# Third interim pass (2026-11-05 03:00): the second pass (interim_analysis2.sh) on the longer trajectories, plus the
# time-resolved FEP convergence diagnostics (a27_fep_convergence.py) on every leg, in parallel.  Same output
# directory ($SCRATCH/tbcb_atomistic/interim/<stamp>/, never read by a21) and the same tgz as the second pass.
# Usage: sbatch -p normal -c 5 --mem=20G -t 1:00:00 --wrap "bash -l scripts/interim_analysis3.sh <stamp> <window_ns>"
module load chemistry
module load py-openmm/8.1.1_py312
module load math
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=1
cd $HOME/tbcb_atomistic
STAMP=$1; W=${2:-110}
O=$SCRATCH/tbcb_atomistic/interim/$STAMP; mkdir -p $O
python3 -W ignore scripts/a27_fep_convergence.py $SCRATCH/tbcb_atomistic/runs/fep --out $O/fep_convergence.json > $O/a27.log 2>&1 &
bash -l scripts/interim_analysis2.sh $STAMP $W
wait
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM3_DONE $O
