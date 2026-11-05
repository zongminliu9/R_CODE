#!/bin/bash
# Marlowe: analysis stages chained (afterany) to the wave-2 compute chains; CPU-only preempt jobs.
export PATH=/cm/shared/apps/slurm/current/bin:$PATH SLURM_CONF=/cm/shared/apps/slurm/var/etc/slurm/slurm.conf
ROOT=/scratch/m000251/zongminl/tbcb_mech/atomistic; cd $ROOT; mkdir -p jobs/pipeline logs
ids() { grep -E "$1" jobs/w2/submitted.txt | cut -d' ' -f2- | tr ' ' '\n' | grep -E '^[0-9]+$' | paste -sd: -; }
sub() { local st=$1 pat=$2; local d=$(ids "$pat"); local dep=""; [ -n "$d" ] && dep="--dependency=afterany:$d"
  sbatch --parsable -p preempt -A marlowe-m000251 -c 4 --mem=16G -t 4:00:00 --exclude=n18,n28 --mail-type=FAIL --mail-user=zongminl@stanford.edu \
    -J ana_$st -o $ROOT/logs/ana_${st}_%j.out $dep $ROOT/scripts/pipeline_stage_marlowe.sh $st; }
f=$(sub fep '^fep_'); r=$(sub rest2 '^rest2_'); m=$(sub md310 '^md310_'); h=$(sub md450 '^md450_')
all=$(cut -d' ' -f2- jobs/w2/submitted.txt | tr ' ' '\n' | grep -E '^[0-9]+$' | paste -sd: -)
rep=$(sbatch --parsable -p preempt -A marlowe-m000251 -c 2 --mem=8G -t 2:00:00 --exclude=n18,n28 --mail-type=END,FAIL --mail-user=zongminl@stanford.edu \
  -J report_tbcb_marlowe -o $ROOT/logs/report_%j.out --dependency=afterany:$all:$f:$r:$m:$h $ROOT/scripts/pipeline_stage_marlowe.sh report)
echo "{\"fep\": \"$f\", \"rest2\": \"$r\", \"md310\": \"$m\", \"md450\": \"$h\", \"report\": \"$rep\"}" | tee jobs/pipeline/pipeline.json
