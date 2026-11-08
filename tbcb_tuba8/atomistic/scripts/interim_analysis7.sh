#!/bin/bash
# Seventh interim pass (2026-11-07): the MD finished after the sixth pass (all four alpha-beta heterodimers, the second
# replicas of the AF3 TBCE complexes, the experimental TBCE complexes threaded on 9M1L) and the final convergence
# diagnostics of every FEP leg after the extensions.
#  - a27 on every FEP leg (folding legs after the S114I extensions, tail legs after the nine extensions).
#  - a26 loop 32-45 CVs and a25 interfaces of the newly finished heterodimer and complex replicas.
# Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/, never read by a21) and a tgz copy.
# Usage: sbatch -p normal -c 4 --mem=32G -t 2:00:00 -o <log> --wrap "bash -l scripts/interim_analysis7.sh <stamp>"
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
(
  $PY scripts/a09_fep_analysis.py $RUNS/fep --out $O/fep_interim.json > $O/a09.log 2>&1
  $PY scripts/a27_fep_convergence.py $RUNS/fep --out $O/fep_convergence.json > $O/a27.log 2>&1
) &
(
  for s in h_9M1M_Tuba1a/r0 h_9M1M_Tuba1a/r1 h_9M1M_Tuba1b/r0 h_9M1M_Tuba1b/r1 h_9M1M_Tuba4a/r0 h_9M1M_Tuba4a/r1 \
           h_9M1M_Tuba8/r0 h_9M1M_Tuba8/r1 c_tbce_a1a_m0/r1 c_tbce_a8_m0/r1 x_9M1L_tbce_a1a/r0 x_9M1L_tbce_a1a/r1 \
           x_9M1L_tbce_a8/r0 x_9M1L_tbce_a8/r1; do
    n=$(echo $s | tr '/' '_')
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_$n.json --stride 2 --dt_ps 20 > $O/a26_md310_$n.log 2>&1
  done
) &
(
  for s in h_9M1M_Tuba1a h_9M1M_Tuba1b h_9M1M_Tuba4a h_9M1M_Tuba8 c_tbce_a1a_m0 c_tbce_a8_m0 x_9M1L_tbce_a1a x_9M1L_tbce_a8; do
    $PY scripts/a25_interface_md.py $RUNS/md310/$s A B $O/iface_$s.json > $O/a25_$s.log 2>&1
  done
) &
wait
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz $STAMP )
echo INTERIM7_DONE $O
