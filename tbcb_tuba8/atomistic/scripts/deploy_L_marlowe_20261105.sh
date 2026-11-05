#!/bin/bash
# Marlowe counterpart of deploy_L_sherlock_20261105.sh: swap the all-L builds from build_L/ into build/, move runs made
# on the defective builds out of runs/, resubmit the two cancelled MD chains (two links each, as in chains_w2.txt),
# release the held jobs, and point the pending analysis/report jobs at the compute jobs that are still queued.
# Usage: bash scripts/deploy_L_marlowe_20261105.sh
set -u
export PATH=/cm/shared/apps/slurm/current/bin:$PATH SLURM_CONF=/cm/shared/apps/slurm/var/etc/slurm/slurm.conf
M=/scratch/m000251/$USER; A=$M/tbcb_mech/atomistic; cd $A
PY=$M/mm/envs/amb/bin/python
NEW="a8 a8_S114I a1a_loop8 a8_A40K a8_I43G a1a_loop8_I114S a8_dom h_9M1M_Tuba8 x_9M1L_tbce_a1a x_9M1L_tbce_a8 c_tbce_a8_m0 t_a1aYF_K188A_cplx"
SUP=$A/superseded_20261105; mkdir -p $SUP/build $SUP/runs/md310 $SUP/runs/rest2 $SUP/runs/fep
for s in $NEW; do [ -s build_L/$s/minimized.cif.gz ] || { echo "missing build_L/$s"; exit 1; }; done
$PY -I scripts/a29_chirality_check.py --strict $(for s in $NEW; do echo build_L/$s/minimized.cif.gz; done) > build_L/a29_strict_deploy.txt
rc=$?; cat build_L/a29_strict_deploy.txt; [ $rc -eq 0 ] || { echo "strict chirality check failed; nothing changed"; exit 1; }
for s in $NEW; do
  [ -d build/$s ] && mv build/$s $SUP/build/$s
  mv build_L/$s build/$s && echo "build/$s replaced"
done
for d in md310/a8 md310/a8_A40K md310/a8_S114I md310/a1a_loop8 rest2/a8_I43G fep/a1a_loop8_114SER_folded; do
  [ -d runs/$d ] && mv runs/$d $SUP/runs/$d && echo "runs/$d moved to superseded"
done
cp -p jobs/w2/submitted.txt $SUP/submitted.txt.before
for j in md310_a8_S114I_r0 md310_a1a_loop8_r0; do
  j1=$(sbatch --parsable jobs/w2/$j.sbatch); j2=$(sbatch --parsable --dependency=afterany:$j1 jobs/w2/$j.sbatch)
  echo "$j -> $j1 $j2"
  grep -v "^$j " jobs/w2/submitted.txt > jobs/w2/submitted.tmp; echo "$j $j1 $j2" >> jobs/w2/submitted.tmp; mv jobs/w2/submitted.tmp jobs/w2/submitted.txt
done
HELD=556001,556002,556003,556004,556013,556014,556023,556024,556025,556047,556048,556051,556052,556055,556056
scontrol release $HELD && echo "released $HELD"
# analysis / report dependencies: every compute job of this campaign that is still queued or running
Q=$(squeue -u $USER -h -o "%i %j")          # one scheduler call, after the resubmissions above
pick() { echo "$Q" | awk -v p="$1" '$2 ~ p {print $1}' | sort -u | paste -sd: -; }
f=$(pick '^fep_'); r=$(pick '^rest2_'); m=$(pick '^md310_'); h=$(pick '^md450_')
P=$(cat jobs/pipeline/pipeline.json)
get() { echo "$P" | sed -E "s/.*\"$1\": \"([0-9]+)\".*/\1/"; }
for st in fep:$f rest2:$r md310:$m md450:$h; do
  k=${st%%:*}; d=${st#*:}; jid=$(get $k)
  [ -n "$d" ] && scontrol update JobId=$jid Dependency=afterany:$d && echo "ana_$k ($jid) -> afterany:$d"
done
all=$(pick '^(md310_|md450_|rest2_|fep_)')
rep=$(get report); scontrol update JobId=$rep Dependency=afterany:$all:$(get fep):$(get rest2):$(get md310):$(get md450) && echo "report ($rep) updated"
