#!/bin/bash
# Eighth (final) interim pass, chained after the last tail analysis (the rerun of TBCE-TUBA1A tail complex r2):
#  - a27 convergence diagnostics of every FEP leg in their final state (a09 is in results/ from ana_tail);
#  - loop 32-45 CVs of the experimental TBCE complexes threaded on 9M1L with the right chain: in their solute.pdb
#    TBCE is chain A and alpha is chain B (pass 7 read chain A, i.e. TBCE residues 32-45; those four rows are void).
# Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/) and a tgz copy.
# Usage: sbatch -p normal -c 2 --mem=16G -t 1:30:00 --dependency=afterany:<ana_tail> -o <log> \
#          --wrap "bash -l scripts/interim_analysis8.sh <stamp>"
module load chemistry
module load py-openmm/8.1.1_py312
module load math
module load py-scikit-learn/1.5.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
cd $HOME/tbcb_atomistic
RUNS=$SCRATCH/tbcb_atomistic/runs
STAMP=$1
O=$SCRATCH/tbcb_atomistic/interim/$STAMP; mkdir -p $O
PY="python3 -W ignore"
$PY scripts/a27_fep_convergence.py $RUNS/fep --out $O/fep_convergence.json > $O/a27.log 2>&1 &
(
  for s in x_9M1L_tbce_a1a/r0 x_9M1L_tbce_a1a/r1 x_9M1L_tbce_a8/r0 x_9M1L_tbce_a8/r1; do
    n=$(echo $s | tr '/' '_')
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_$n.json --chain B --stride 2 --dt_ps 20 > $O/a26_md310_$n.log 2>&1
  done
) &
wait
cp -p $SCRATCH/tbcb_atomistic/results/fep_results.json $O/fep_interim.json
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM8_DONE $O
