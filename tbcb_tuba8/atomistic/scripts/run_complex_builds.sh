#!/bin/bash
cd /scratch/m000251/zongminl/tbcb_mech/atomistic
PY=/scratch/m000251/zongminl/mm/envs/amb/bin/python
export OPENMM_CPU_THREADS=4 OMP_NUM_THREADS=4
A=inputs/af3
b() { sid=$1; shift; [ -s build/$sid/build.json ] || [ -s build/$sid/vacuum_forces.json -a "$VAC" = 1 ] || $PY scripts/a07_build_complex.py $sid "$@" > logs/cbuild_$sid.log 2>&1; echo "$(date) $sid" >> logs/run_complex_builds.progress; }
VAC=0
b c_tbcb_a8_m0 af3 $A/fold_tbcb_tuba8/fold_tbcb_tuba8_model_0.cif --alpha B --partners A
b c_tbce_a1a_m0 af3 $A/fold_tbce_tuba1a/fold_tbce_tuba1a_model_0.cif --alpha A --partners B
b c_tbce_a8_m0 af3 $A/fold_tbce_tuba8/fold_tbce_tuba8_model_0.cif --alpha A --partners B
b x_9M1L_tbce_a1a exp inputs/9M1L.cif.gz --alpha a --partners E --thread Tuba1a
b x_9M1L_tbce_a8 exp inputs/9M1L.cif.gz --alpha a --partners E --thread Tuba8
VAC=1
for k in 1 2 3 4; do
  b v_tbcb_a1a_m$k af3 $A/fold_tbcb_tuba1a/fold_tbcb_tuba1a_model_$k.cif --alpha B --partners A --vacuum-only
  b v_tbcb_a8_m$k af3 $A/fold_tbcb_tuba8/fold_tbcb_tuba8_model_$k.cif --alpha B --partners A --vacuum-only
  b v_tbce_a1a_m$k af3 $A/fold_tbce_tuba1a/fold_tbce_tuba1a_model_$k.cif --alpha A --partners B --vacuum-only
  b v_tbce_a8_m$k af3 $A/fold_tbce_tuba8/fold_tbce_tuba8_model_$k.cif --alpha A --partners B --vacuum-only
done
echo ALLDONE >> logs/run_complex_builds.progress
