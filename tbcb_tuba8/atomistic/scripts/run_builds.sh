#!/bin/bash
# build + minimise every monomer construct on the Marlowe login node, 4 threads, sequential
cd /scratch/m000251/zongminl/tbcb_mech/atomistic
PY=/scratch/m000251/zongminl/mm/envs/amb/bin/python
export OPENMM_CPU_THREADS=4 OMP_NUM_THREADS=4
while IFS=$'\t' read -r sid spec; do
  [ -z "$sid" ] && continue
  if [ ! -s build/$sid/solvated.cif.gz ]; then $PY scripts/a01_build.py $sid "$spec" > logs/build_$sid.log 2>&1; fi
  if [ ! -s build/$sid/minimize.json ] || ! [ -f build/$sid/minimize.flag ]; then
     $PY scripts/a02_minimize.py build/$sid > logs/min_$sid.log 2>&1 && echo "flexible" >> build/$sid/minimize.flag
  fi
  echo "$(date) $sid done" >> logs/run_builds.progress
done < scripts/systems_monomer.tsv
echo ALLDONE >> logs/run_builds.progress
