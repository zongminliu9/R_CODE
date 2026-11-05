#!/bin/bash
# Fourth interim pass (2026-11-05, after md310_a8_r0 reaches 200 ns): FEP progress (a09, a27), REST2 progress (a10),
# loop CVs of the four free monomers, and the TUBA8 monomer residue level before and after the loop unfolds
# (10-140 ns intermediate, 140-200 ns unfolded), computed with a13_residue_analysis_next.py, which adds the
# periodic-image check for GTP / Mg / Mg-shell waters (this run is also its test before it replaces a13).
# It also repairs the five builds with wrong chirality (a30_fix_chirality.py: four mirror-image complexes from the
# a07 principal-axes bug are reflected back, t_a1aYF_K188A_cplx is re-minimised), re-checks every build with
# a29_chirality_check.py, and runs a26 (now handedness-aware) on the mirror-image c_tbcb_a8_m0 r0 trajectory and
# on an L control.  The ten jobs that use the five builds are held until this has been checked.
# Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/, never read by a21) and a tgz copy.
# Usage: sbatch -p normal -c 6 --mem=32G -t 1:45:00 --dependency=afterany:<md310_a8 job> \
#          --wrap "bash -l scripts/interim_analysis4.sh <stamp>"
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
  $PY scripts/a30_fix_chirality.py $BUILD > $O/a30.log 2>&1
  $PY scripts/a29_chirality_check.py $BUILD/*/minimized.cif.gz > $O/a29_builds_after_fix.txt 2>&1
  echo "a29 exit $?" >> $O/a29_builds_after_fix.txt
  # release the ten jobs held for the repair only if the repair finished and every build is L
  HELD=51130003,51130004,51131394,51130007,51130008,51131979,51131981,51131984,51131986,51130634
  if grep -q CHIRALITY_FIX_DONE $O/a30.log && grep -qx "a29 exit 0" $O/a29_builds_after_fix.txt; then
    scontrol release $HELD && echo "released $HELD" > $O/release.txt
  else
    echo "repair or check failed; jobs stay held" > $O/release.txt
  fi
) &
$PY scripts/a13_residue_analysis_next.py $RUNS/md310/a8 $BUILD/a8 $O/residue_a8_10-140ns.json --skip_ns 10 --max_ns 140 --stride 5 > $O/a13_a8_10-140.log 2>&1 &
$PY scripts/a13_residue_analysis_next.py $RUNS/md310/a8 $BUILD/a8 $O/residue_a8_140-200ns.json --skip_ns 140 --max_ns 200 --stride 5 > $O/a13_a8_140-200.log 2>&1 &
(
  $PY scripts/a09_fep_analysis.py $RUNS/fep --out $O/fep_interim.json > $O/a09.log 2>&1
  $PY scripts/a27_fep_convergence.py $RUNS/fep --out $O/fep_convergence.json > $O/a27.log 2>&1
) &
(
  for d in $RUNS/rest2/*/r0; do sid=$(basename $(dirname $d))
    $PY scripts/a10_rest2_analysis.py $d --thresholds $SUM/rest2_state_thresholds.json > $O/a10_$sid.log 2>&1 && cp $d/rest2_results.json $O/rest2_$sid.json
  done
  for s in a1a/r0 a1b/r0 a4a/r0 a8/r0 c_tbcb_a8_m0/r0 h_9M1M_Tuba1a/r1; do
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_${s%/*}.json --stride 2 --dt_ps 20 > $O/a26_md310_${s%/*}.log 2>&1
  done
  for sid in a1a a8 a1a_loop8 a8_loop1A; do
    $PY scripts/a26_loop_cvs.py $RUNS/rest2/$sid/r0 $O/loopcv_rest2s0_$sid.json --glob 'rest2_s0_from*.dcd' --stride 5 --dt_ps 2 > $O/a26_rest2_$sid.log 2>&1
  done
) &
wait
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM4_DONE $O
