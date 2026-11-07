#!/bin/bash
# Sixth interim pass (2026-11-06): everything that finished after the fifth pass, recomputed from the complete run files.
#  - FEP: a09 + a27 on every leg (folding legs: the final I114S legs; tail legs: Y451F TBCB-TUBA1A, F449Y TBCB-TUBA8
#    with and without K188A, F449Y TBCE-TUBA8 complete; two more part-way), and a22 (ddG_bind, Kd ratio, tauRAMD
#    residence times) on a copy of this pass's a09 output with the ramd runs linked in.
#  - md310: loop 32-45 CVs (a26) of the four free monomers (200 ns each; TUBA8 rerun on the all-L build) and of every
#    finished complex / heterodimer replica; interfaces (a25) of the complexes (replicas still running are read as far
#    as they have got and are labelled by their trajectory length).
#  - US-phi43: a17 analyse for the three finished runs (TUBA1A, TUBA1A-G43I, TUBA1A-loop8).
#  - MSM: a11 with the deeptime 0.4.5 ck_test call fixed, written here (not to results/) as its test.
# Output: $SCRATCH/tbcb_atomistic/interim/<stamp>/ (outside results/, never read by a21) and a tgz copy.
# Usage: sbatch -p normal -c 6 --mem=48G -t 2:30:00 -o <log> --wrap "bash -l scripts/interim_analysis6.sh <stamp>"
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
  mkdir -p $O/binding && cp $O/fep_interim.json $O/binding/fep_results.json && ln -sfn $RUNS/ramd $O/binding/ramd
  $PY scripts/a22_binding_analysis.py $O/binding > $O/a22.log 2>&1
  echo "a22 exit $?" >> $O/a22.log
  $PY scripts/a27_fep_convergence.py $RUNS/fep --out $O/fep_convergence.json > $O/a27.log 2>&1
) &
(
  for s in a1a/r0 a1b/r0 a4a/r0 a8/r0 h_9M1M_Tuba1a/r1 h_9M1M_Tuba4a/r0 c_tbcb_a1a_m0/r0 c_tbcb_a1a_m0/r1 c_tbcb_a1a_m0/r2 \
           c_tbcb_a8_m0/r0 c_tbcb_a8_m0/r1 c_tbcb_a8_m0/r2 c_tbce_a1a_m0/r0 c_tbce_a8_m0/r0; do
    n=$(echo $s | tr '/' '_')
    $PY scripts/a26_loop_cvs.py $RUNS/md310/$s $O/loopcv_md310_$n.json --stride 2 --dt_ps 20 > $O/a26_md310_$n.log 2>&1
  done
) &
(
  $PY scripts/a25_interface_md.py $RUNS/md310/c_tbcb_a1a_m0 A B $O/iface_c_tbcb_a1a_m0.json \
      --focus B:187,B:188,B:192,B:203,B:204,B:223,A:448,A:449,A:450,A:451 > $O/a25_c_tbcb_a1a_m0.log 2>&1
  $PY scripts/a25_interface_md.py $RUNS/md310/c_tbcb_a8_m0 A B $O/iface_c_tbcb_a8_m0.json \
      --focus B:187,B:188,B:192,B:203,B:204,B:223,A:446,A:447,A:448,A:449 > $O/a25_c_tbcb_a8_m0.log 2>&1
  for s in c_tbce_a1a_m0 c_tbce_a8_m0 h_9M1M_Tuba4a h_9M1M_Tuba1a; do
    $PY scripts/a25_interface_md.py $RUNS/md310/$s A B $O/iface_$s.json > $O/a25_$s.log 2>&1
  done
) &
(
  for s in a1a a1a_G43I a1a_loop8; do
    $PY scripts/a17_us_phi.py analyse $RUNS/usphi43/$s/r0 > $O/a17_$s.log 2>&1 && cp $RUNS/usphi43/$s/r0/usphi_results.json $O/usphi_$s.json
  done
  $PY scripts/a11_msm.py $RUNS/md310 $O/msm --systems a1a,a1b,a4a,a8 > $O/a11.log 2>&1
  echo "a11 exit $?" >> $O/a11.log
) &
wait
ls -la $O
( cd $SCRATCH/tbcb_atomistic/interim && tar czf $HOME/tbcb_atomistic/results_summary/interim_$STAMP.tgz --exclude='binding/ramd' $STAMP )
echo INTERIM6_DONE $O
