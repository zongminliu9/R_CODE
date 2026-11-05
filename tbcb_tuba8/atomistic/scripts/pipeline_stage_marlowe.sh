#!/bin/bash
# Analysis stage run inside a Slurm job that depends (afterany) on the compute jobs of its group.
# Usage: pipeline_stage.sh <stage>   stages: fep rest2 md310 md450 pmf usphi tail mmgbsa report
# Every sub-step is attempted even if another fails; results/stage_<stage>.json records which succeeded.
STAGE=$1




PYAMB=/scratch/m000251/zongminl/mm/envs/amb/bin/python
export OPENMM_CPU_THREADS=${SLURM_CPUS_PER_TASK:-4}
ROOT=/scratch/m000251/zongminl/tbcb_mech/atomistic; cd $ROOT
RUNS=$ROOT/runs
RES=$ROOT/results
SUM=$ROOT/results_summary
mkdir -p $RES $SUM
PY="$PYAMB -W ignore"
ok=(); bad=()
step() { local name=$1; shift; echo "== $name: $*"; if "$@"; then ok+=("$name"); else bad+=("$name"); fi; }
case $STAGE in
  fep)
    step a09 $PY scripts/a09_fep_analysis.py $RUNS/fep --out $RES/fep_results.json ;;
  rest2)
    for d in $RUNS/rest2/*/r0; do sid=$(basename $(dirname $d)); [ -f $d/rest2_ck.npz ] || continue
      mkdir -p $RES/rest2/$sid
      step a10_$sid $PY scripts/a10_rest2_analysis.py $d --thresholds $SUM/rest2_state_thresholds.json && cp $d/rest2_results.json $RES/rest2/$sid/
      step a23_$sid $PY scripts/a23_representatives.py $d $RES/representatives/$sid
    done ;;
  md310)
    for d in $RUNS/md310/a*; do sid=$(basename $d); mkdir -p $RES/residue
      step a13_$sid $PY scripts/a13_residue_analysis.py $d $ROOT/build/$sid $RES/residue/$sid.json
    done
    sys=$(ls -d $RUNS/md310/a* 2>/dev/null | xargs -n1 basename | grep -v _dom | paste -sd, -)
    [ -n "$sys" ] && step a11_msm $PY scripts/a11_msm.py $RUNS/md310 $RES/msm --systems $sys ;;
  md450)
    for d in $RUNS/md450/*; do sid=$(basename $d); mkdir -p $RES/unfold
      step a18_$sid $PY scripts/a18_unfold_analysis.py $d $ROOT/build/$sid $RES/unfold/$sid.json
    done ;;
  pmf)
    for kind in pmf pmfdom; do for d in $RUNS/$kind/*/r0; do sid=$(basename $(dirname $d)); [ -f $d/pmf_ck.npz ] || continue
      mkdir -p $RES/$kind/$sid
      step a15_$sid $PY scripts/a15_pmf_tail.py analyse $d && cp $d/pmf_results.json $RES/$kind/$sid/
    done; done ;;
  usphi)
    for d in $RUNS/usphi43/*/r0; do sid=$(basename $(dirname $d)); [ -f $d/usphi_ck.npz ] || continue
      mkdir -p $RES/usphi43/$sid
      step a17_$sid $PY scripts/a17_us_phi.py analyse $d && cp $d/usphi_results.json $RES/usphi43/$sid/
    done ;;
  tail)
    step a09 $PY scripts/a09_fep_analysis.py $RUNS/fep --out $RES/fep_results.json
    [ -e $RES/ramd ] || ln -s $RUNS/ramd $RES/ramd
    step a22 $PY scripts/a22_binding_analysis.py $RES ;;
  mmgbsa)
    for p in $ROOT/mmpbsa/*; do sid=$(basename $p); [ -d $RUNS/md310/$sid ] || continue
      mkdir -p $RES/mmpbsa
      step a12b_$sid $PY scripts/a12b_mmgbsa_openmm.py $p $RUNS/md310/$sid $RES/mmpbsa/$sid.json --platform ${MMPLAT:-CUDA}
    done ;;
  report)
    step a09 $PY scripts/a09_fep_analysis.py $RUNS/fep --out $RES/fep_results.json
    [ -e $RES/ramd ] || ln -s $RUNS/ramd $RES/ramd
    step a22 $PY scripts/a22_binding_analysis.py $RES
    [ -f $SUM/minimization_forces.json ] && cp $SUM/minimization_forces.json $RES/
    step a21 $PY scripts/a21_compile.py $RES $RES/master
    step a24 $PY scripts/a24_report.py $RES $RES/report_zh.md ;;
esac
# small, purge-proof copy of every summary JSON / CSV / report under $HOME
( cd $RES && find . -maxdepth 3 \( -name "*.json" -o -name "*.csv" -o -name "*.md" \) -size -20M -print0 | tar czf $SUM/results_summary_$STAGE.tgz --null -T - ) 2>/dev/null
cp -f $RES/master.csv $RES/master.json $RES/report_zh.md $SUM/ 2>/dev/null
jl() { [ $# -eq 0 ] && return; printf '"%s",' "$@" | sed 's/,$//'; }
printf '{"stage": "%s", "time": "%s", "job": "%s", "ok": [%s], "failed": [%s]}\n' "$STAGE" "$(date)" "$SLURM_JOB_ID" \
  "$(jl "${ok[@]}")" "$(jl "${bad[@]}")" > $RES/stage_$STAGE.json
cp $RES/stage_$STAGE.json $SUM/
cat $RES/stage_$STAGE.json
[ ${#bad[@]} -eq 0 ]
