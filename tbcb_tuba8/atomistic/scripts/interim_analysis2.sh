#!/bin/bash
# Second interim pass (2026-11-05): FEP MBAR including the completed I114S legs, REST2 with the calibrated loop states
# (a10 --thresholds), residue analysis with the SASA unit fixed on a common window, interface and loop CVs on the
# longer trajectories.  Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/, never read by a21).
# Chain letters are those of solute.pdb: alpha = A in every complex, TBCB / TBCE / beta = B.
# Usage: sbatch -p normal -c 4 --mem=16G -t 1:00:00 --wrap "bash -l scripts/interim_analysis2.sh <stamp> <window_ns>"
module load chemistry
module load py-openmm/8.1.1_py312
module load math
module load py-scikit-learn/1.5.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
cd $HOME/tbcb_atomistic
RUNS=$SCRATCH/tbcb_atomistic/runs; BUILD=$SCRATCH/tbcb_atomistic/build; SUM=$HOME/tbcb_atomistic/results_summary
STAMP=$1; W=${2:-85}
O=$SCRATCH/tbcb_atomistic/interim/$STAMP; mkdir -p $O
PY="python3 -W ignore"
(
  $PY scripts/a09_fep_analysis.py $RUNS/fep --out $O/fep_interim.json > $O/a09.log 2>&1
  for d in $RUNS/rest2/*/r0; do sid=$(basename $(dirname $d))
    $PY scripts/a10_rest2_analysis.py $d --thresholds $SUM/rest2_state_thresholds.json > $O/a10_$sid.log 2>&1 && cp $d/rest2_results.json $O/rest2_$sid.json
  done
  for s in a1a/r0 a1b/r0 a4a/r0 a8/r0 h_9M1M_Tuba1a/r1 c_tbcb_a8_m0/r0 c_tbce_a8_m0/r0; do
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_${s%/*}.json --stride 2 --dt_ps 20 > $O/a26_md310_${s%/*}.log 2>&1
  done
  for sid in a1a a8 a1a_loop8 a8_loop1A; do
    $PY scripts/a26_loop_cvs.py $RUNS/rest2/$sid/r0 $O/loopcv_rest2s0_$sid.json --glob 'rest2_s0_from*.dcd' --stride 5 --dt_ps 2 > $O/a26_rest2_$sid.log 2>&1
  done
) &
for sid in a1a a1b a4a; do
  $PY scripts/a13_residue_analysis.py $RUNS/md310/$sid $BUILD/$sid $O/residue_${sid}_10-${W}ns.json --max_ns $W --stride 5 > $O/a13_${sid}.log 2>&1 &
done
wait
$PY scripts/a13_residue_analysis.py $RUNS/md310/a8 $BUILD/a8 $O/residue_a8_10-${W}ns.json --max_ns $W --stride 5 > $O/a13_a8.log 2>&1 &
$PY scripts/a25_interface_md.py $RUNS/md310/c_tbcb_a8_m0 A B $O/iface_c_tbcb_a8_m0.json --focus B:187,B:188,B:192,B:203,B:204,B:223,A:446,A:447,A:448,A:449 > $O/a25_c_tbcb_a8_m0.log 2>&1 &
$PY scripts/a25_interface_md.py $RUNS/md310/c_tbce_a8_m0 A B $O/iface_c_tbce_a8_m0.json > $O/a25_c_tbce_a8_m0.log 2>&1 &
$PY scripts/a25_interface_md.py $RUNS/md310/h_9M1M_Tuba1a A B $O/iface_h_9M1M_Tuba1a.json > $O/a25_h_9M1M_Tuba1a.log 2>&1 &
wait
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM_DONE $O
