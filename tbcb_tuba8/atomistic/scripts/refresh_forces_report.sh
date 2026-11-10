#!/bin/bash
# The minimisation / residual-force table in results/ was written on 2026-11-04, before twelve builds were rebuilt
# all-L on 2026-11-05, so the final report listed forces of the superseded builds.  Rebuild the table from the current
# build/ directory, then rerun the report stage (a09, a22, a21, a24).
# Usage: sbatch -p normal -c 2 --mem=8G -t 1:00:00 -J report_tbcb -o <log> --wrap "bash -l scripts/refresh_forces_report.sh"
module load chemistry
module load py-openmm/8.1.1_py312
export PYTHONPATH=$PYTHONPATH:$HOME/fep_pylib
cd $HOME/tbcb_atomistic
RES=$SCRATCH/tbcb_atomistic/results
mkdir -p $SCRATCH/tbcb_atomistic/superseded_20261107/results
[ -f $SCRATCH/tbcb_atomistic/superseded_20261107/results/minimization_forces_20261104.json ] || cp -p $RES/minimization_forces.json $SCRATCH/tbcb_atomistic/superseded_20261107/results/minimization_forces_20261104.json
python3 -W ignore scripts/a20_force_table.py $SCRATCH/tbcb_atomistic/build $RES/minimization_forces > $SCRATCH/tbcb_atomistic/logs/a20_refresh.log 2>&1
echo "a20 exit $?"
SUM=$HOME/tbcb_atomistic/results_summary
for f in minimization_forces.json minimization_forces.csv; do
  [ -f $SUM/$f ] && [ ! -f $SCRATCH/tbcb_atomistic/superseded_20261107/results/summary_$f ] && cp -p $SUM/$f $SCRATCH/tbcb_atomistic/superseded_20261107/results/summary_$f
  cp -p $RES/$f $SUM/$f
done
bash scripts/pipeline_stage.sh report
