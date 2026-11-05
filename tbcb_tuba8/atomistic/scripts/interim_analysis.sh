#!/bin/bash
# Interim (mid-run) analysis of what the running Sherlock jobs have written so far; nothing here feeds the final
# pipeline.  Run directories are only read, except a10, which writes rest2_results.json next to its checkpoint
# (the final ana_rest2 stage overwrites it).  Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/,
# so a21_compile never picks it up) and a tgz copy under $HOME/tbcb_atomistic/results_summary/.
# Chain letters below are those of the solute.pdb files (OpenMM renamed chains in order): alpha = A in every
# complex; TBCB / TBCE / beta = B.
# Usage: sbatch -p dev -c 4 --mem=16G -t 1:00:00 scripts/interim_analysis.sh [stamp] [w1_ns] [w2_ns]
module load chemistry
module load py-openmm/8.1.1_py312
module load math
module load py-scikit-learn/1.5.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
cd $HOME/tbcb_atomistic
RUNS=$SCRATCH/tbcb_atomistic/runs; BUILD=$SCRATCH/tbcb_atomistic/build
STAMP=${1:-$(date +%Y%m%d_%H%M)}; W1=${2:-24}; W2=${3:-41}
O=$SCRATCH/tbcb_atomistic/interim/$STAMP; mkdir -p $O
PY="python3 -W ignore"
# loop 32-45 CVs recomputed from the 310 K trajectories (calibration of the REST2 state definitions)
(
  for s in a1a/r0 a1b/r0 a4a/r0 a8/r0 h_9M1M_Tuba1a/r1 c_tbcb_a8_m0/r0 c_tbce_a8_m0/r0; do
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_${s%/*}.json --stride 2 --dt_ps 20 > $O/a26_md310_${s%/*}.log 2>&1
  done
  for sid in a1a a8 a1a_loop8 a8_loop1A; do
    $PY scripts/a26_loop_cvs.py $RUNS/rest2/$sid/r0 $O/loopcv_rest2s0_$sid.json --glob 'rest2_s0_from*.dcd' --stride 5 --dt_ps 2 > $O/a26_rest2_$sid.log 2>&1
  done
) &
(
  $PY scripts/a09_fep_analysis.py $RUNS/fep --out $O/fep_interim.json > $O/a09.log 2>&1
  for d in $RUNS/rest2/*/r0; do sid=$(basename $(dirname $d))
    $PY scripts/a10_rest2_analysis.py $d > $O/a10_$sid.log 2>&1 && cp $d/rest2_results.json $O/rest2_$sid.json
  done
) &
for sid in a1a a1b a4a; do
  $PY scripts/a13_residue_analysis.py $RUNS/md310/$sid $BUILD/$sid $O/residue_${sid}_10-${W1}ns.json --max_ns $W1 --stride 2 > $O/a13_${sid}_$W1.log 2>&1 &
done
wait
$PY scripts/a13_residue_analysis.py $RUNS/md310/a8 $BUILD/a8 $O/residue_a8_10-${W1}ns.json --max_ns $W1 --stride 2 > $O/a13_a8_$W1.log 2>&1 &
$PY scripts/a25_interface_md.py $RUNS/md310/c_tbcb_a8_m0 A B $O/iface_c_tbcb_a8_m0.json --focus B:187,B:188,B:203,B:204,A:447,A:448,A:449 > $O/a25_c_tbcb_a8_m0.log 2>&1 &
$PY scripts/a25_interface_md.py $RUNS/md310/c_tbce_a8_m0 A B $O/iface_c_tbce_a8_m0.json --focus B:31,B:32,B:48,B:49,B:68,B:69,A:447,A:448,A:449 > $O/a25_c_tbce_a8_m0.log 2>&1 &
$PY scripts/a25_interface_md.py $RUNS/md310/h_9M1M_Tuba1a A B $O/iface_h_9M1M_Tuba1a.json > $O/a25_h_9M1M_Tuba1a.log 2>&1 &
wait
for sid in a1a a1b a8; do
  $PY scripts/a13_residue_analysis.py $RUNS/md310/$sid $BUILD/$sid $O/residue_${sid}_10-${W2}ns.json --max_ns $W2 --stride 2 > $O/a13_${sid}_$W2.log 2>&1 &
done
wait
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM_DONE $O
