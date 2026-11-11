#!/bin/bash
# Marlowe CPU job (2026-11-10): analysis of the Marlowe runs that finished after 11-06 (about 5 CPU-h on preempt):
#  - FEP a09 + a27 on every Marlowe leg (I110V, A270V, A294S, S114I in the TUBA8-with-TUBA1A-loop context complete;
#    I114S in the TUBA1A-with-TUBA8-loop context 2 of 3 complete);
#  - loop 32-45 CVs (a26) and residue level (a13) of the finished 200 ns runs (TUBA1A r1/r2, K40A, both loop swaps,
#    TUBA1A-I114S at 197.7 ns);
#  - 450 K unfolding (a18) of TUBA1A-I114S, TUBA1B, TUBA4A, TUBA8-S114I (3 x 20 ns each);
#  - REST2 (a10, calibrated states) of TUBA1A-G43I and TUBA1B (30 ns per replica).
# Output: interim/<stamp>/ (never results/).
# Usage: sbatch -p preempt -A marlowe-m000251 -c 8 --mem=32G -t 2:00:00 -o <log> scripts/interim_marlowe2.sh <stamp>
M=/scratch/m000251/$USER
cd $M/tbcb_mech/atomistic
PY="$M/mm/envs/amb/bin/python -W ignore"
export OMP_NUM_THREADS=1 OPENMM_CPU_THREADS=1
O=interim/$1; mkdir -p $O
( $PY scripts/a09_fep_analysis.py runs/fep --out $O/fep_interim.json > $O/a09.log 2>&1
  $PY scripts/a27_fep_convergence.py runs/fep --out $O/fep_convergence.json > $O/a27.log 2>&1 ) &
( for s in a1a/r1 a1a/r2 a1a_K40A/r0 a1a_loop8/r0 a8_loop1A/r0 a1a_I114S/r0; do n=$(echo $s | tr '/' '_')
    $PY scripts/a26_loop_cvs.py runs/md310/$s $O/loopcv_md310_$n.json --stride 2 --dt_ps 20 > $O/a26_$n.log 2>&1; done
  for s in a1a_I114S a1b a4a a8_S114I; do
    $PY scripts/a18_unfold_analysis.py runs/md450/$s build/$s $O/unfold_md450_$s.json > $O/a18_$s.log 2>&1; done
  for s in a1a_G43I a1b; do
    $PY scripts/a10_rest2_analysis.py runs/rest2/$s/r0 > $O/a10_$s.log 2>&1 && cp runs/rest2/$s/r0/rest2_results.json $O/rest2_$s.json; done ) &
for s in a1a a1a_K40A a1a_loop8 a8_loop1A a1a_I114S; do
  $PY scripts/a13_residue_analysis.py runs/md310/$s build/$s $O/residue_${s}_10-200ns.json --skip_ns 10 --max_ns 200 --stride 5 > $O/a13_$s.log 2>&1 &
done
wait
ls -la $O
echo INTERIM_MARLOWE2_DONE
