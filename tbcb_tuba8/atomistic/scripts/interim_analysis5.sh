#!/bin/bash
# Fifth interim pass (2026-11-05 07:00): everything that finished after the fourth pass, recomputed from the complete
# run files rather than from logs.
#  - FEP: a09 + a27 on every leg.  The three Y451F legs of the TBCB-TUBA1A tail complex finished at 02:47 / 05:02 /
#    05:22; the three I114S extension legs are part-way, so their values here are progress values, not final.
#  - md310: loop CVs (a26) of the free monomers and of the two finished complexes; residue level (a13, with the
#    periodic-image check) of TUBA1A / TUBA1B / TUBA4A on a common 10-150 ns window, and of TUBA4A on 10-200 ns
#    (its 200 ns run finished at 05:19).
#  - md450: unfolding analysis (a18) of the three finished 450 K TUBA1A replicas (3 x 20 ns).
#  - REST2 progress (a10), and a strict per-residue chirality check of the solutes written by the runs that started
#    on the rebuilt all-L builds.
# Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/, never read by a21) and a tgz copy.
# Usage: sbatch -p normal -c 6 --mem=40G -t 1:45:00 -o <log> --wrap "bash -l scripts/interim_analysis5.sh <stamp>"
module load chemistry
module load py-openmm/8.1.1_py312
module load math
module load py-scikit-learn/1.5.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
cd $HOME/tbcb_atomistic
RUNS=$SCRATCH/tbcb_atomistic/runs; BUILD=$SCRATCH/tbcb_atomistic/build; SUM=$HOME/tbcb_atomistic/results_summary
STAMP=$1
O=$SCRATCH/tbcb_atomistic/interim/$STAMP; mkdir -p $O
PY="python3 -W ignore"
(
  $PY scripts/a09_fep_analysis.py $RUNS/fep --out $O/fep_interim.json > $O/a09.log 2>&1
  $PY scripts/a27_fep_convergence.py $RUNS/fep --out $O/fep_convergence.json > $O/a27.log 2>&1
) &
$PY scripts/a13_residue_analysis.py $RUNS/md310/a4a $BUILD/a4a $O/residue_a4a_10-200ns.json --skip_ns 10 --max_ns 200 --stride 5 > $O/a13_a4a_10-200.log 2>&1 &
for sid in a1a a1b a4a; do
  $PY scripts/a13_residue_analysis.py $RUNS/md310/$sid $BUILD/$sid $O/residue_${sid}_10-150ns.json --skip_ns 10 --max_ns 150 --stride 5 > $O/a13_${sid}_10-150.log 2>&1 &
done
(
  $PY scripts/a18_unfold_analysis.py $RUNS/md450/a1a $BUILD/a1a $O/unfold_md450_a1a.json > $O/a18_md450_a1a.log 2>&1
  for s in a1a/r0 a1b/r0 a4a/r0 c_tbcb_a8_m0/r0 h_9M1M_Tuba1a/r1; do
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_${s%/*}.json --stride 2 --dt_ps 20 > $O/a26_md310_${s%/*}.log 2>&1
  done
  for d in $RUNS/rest2/*/r0; do sid=$(basename $(dirname $d))
    $PY scripts/a10_rest2_analysis.py $d --thresholds $SUM/rest2_state_thresholds.json > $O/a10_$sid.log 2>&1 && cp $d/rest2_results.json $O/rest2_$sid.json
  done
  NEWSOL=$(ls $RUNS/md450/a8/r*/solute.pdb $RUNS/rest2/a8/r0/solute.pdb $RUNS/md310/a8/r*/solute.pdb $RUNS/md310/c_tbce_a8_m0/r*/solute.pdb 2>/dev/null)
  $PY scripts/a29_chirality_check.py --strict $NEWSOL > $O/a29_new_runs.txt 2>&1
  echo "a29 exit $?" >> $O/a29_new_runs.txt
) &
wait
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM5_DONE $O
