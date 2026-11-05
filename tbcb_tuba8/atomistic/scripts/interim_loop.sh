#!/bin/bash
# Loop 32-45 CVs recomputed from the 310 K trajectories (MD and REST2 state 0) to calibrate the REST2 state
# definitions; same output directory convention as interim_analysis.sh.  Usage: interim_loop.sh <stamp>
module load chemistry
module load py-openmm/8.1.1_py312
module load math
module load py-scikit-learn/1.5.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
export OMP_NUM_THREADS=1
cd $HOME/tbcb_atomistic
RUNS=$SCRATCH/tbcb_atomistic/runs; O=$SCRATCH/tbcb_atomistic/interim/$1; mkdir -p $O
PY="python3 -W ignore"
for s in a1a/r0 a1b/r0 a4a/r0 a8/r0 h_9M1M_Tuba1a/r1 c_tbcb_a8_m0/r0 c_tbce_a8_m0/r0; do
  $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_${s%/*}.json --stride 2 --dt_ps 20 > $O/a26_md310_${s%/*}.log 2>&1 &
done
wait
for sid in a1a a8 a1a_loop8 a8_loop1A; do
  $PY scripts/a26_loop_cvs.py $RUNS/rest2/$sid/r0 $O/loopcv_rest2s0_$sid.json --glob 'rest2_s0_from*.dcd' --stride 5 --dt_ps 2 > $O/a26_rest2_$sid.log 2>&1 &
done
wait
ls $O/loopcv_* ; echo LOOP_DONE
