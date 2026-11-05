#!/bin/bash
# Swap the rebuilt (all-L) builds into $SCRATCH/tbcb_atomistic/build on Sherlock, move the runs made on the defective
# builds out of runs/, resubmit the stopped / finished-but-invalid runs, and release the held jobs.  Light file
# operations and scheduler calls only (login node).  Expects the new builds unpacked in $SCRATCH/tbcb_atomistic/build_L
# and a29 --strict to have passed on them (checked again here; nothing is changed if it fails).
# Usage: bash scripts/deploy_L_sherlock_20261105.sh            (prints what it does)
set -u
A=$SCRATCH/tbcb_atomistic; cd $HOME/tbcb_atomistic
NEW="a8 a8_S114I a1a_loop8 a8_A40K a8_I43G a1a_loop8_I114S a8_dom h_9M1M_Tuba8 x_9M1L_tbce_a1a x_9M1L_tbce_a8 c_tbce_a8_m0 t_a1aYF_K188A_cplx"
SUP=$A/superseded_20261105; mkdir -p $SUP/build $SUP/runs/md310 $SUP/runs/rest2
module load chemistry 2>/dev/null; module load py-openmm/8.1.1_py312 2>/dev/null
for s in $NEW; do [ -s $A/build_L/$s/minimized.cif.gz ] || { echo "missing build_L/$s"; exit 1; }; done
python3 scripts/a29_chirality_check.py --strict $(for s in $NEW; do echo $A/build_L/$s/minimized.cif.gz; done) > $A/build_L/a29_strict_sherlock.txt
rc=$?; cat $A/build_L/a29_strict_sherlock.txt; [ $rc -eq 0 ] || { echo "strict chirality check failed; nothing changed"; exit 1; }
# builds: old -> superseded, new -> build (one rename each, so a job starting now sees either the old or the new build)
for s in $NEW; do
  [ -d $A/build/$s ] && mv $A/build/$s $SUP/build/$s
  mv $A/build_L/$s $A/build/$s && echo "build/$s replaced"
done
# runs made on the defective builds (the S114I FEP legs stay: they use their own run-directory copies)
for d in md310/a8 md310/c_tbce_a8_m0 rest2/a8 rest2/a1a_loop8; do
  [ -d $A/runs/$d ] && mv $A/runs/$d $SUP/runs/$d && echo "runs/$d moved to superseded"
done
# resubmit
for j in w1/md310_a8_r0 w1/rest2_a8_r0 w1/rest2_a1a_loop8_r0 w3/md310_c_tbce_a8_m0_r0; do
  echo "$j -> $(sbatch --parsable jobs/$j.sbatch)"
done
# release exactly the jobs held on 2026-11-05 for the chirality repair (all of them now use L builds)
HELD=51130003,51130004,51131394,51130007,51130008,51131979,51131981,51131984,51131986,51130634,51132836,51132821,51132829,51132967,51131993,51131995,51130012
scontrol release $HELD && echo "released $HELD"
python3 scripts/pipeline_submit.py | tail -12
