#!/bin/bash
# Marlowe CPU job: FEP analysis (a09 + a27) of the legs run on Marlowe (A270V and A294S complete, others part-way),
# written to an interim directory (never to results/).  About 4 CPU-h on the preempt partition.
# Usage: sbatch -p preempt -A marlowe-m000251 -c 8 --mem=16G -t 1:00:00 -o <log> scripts/interim_marlowe_fep.sh <stamp>
M=/scratch/m000251/$USER
cd $M/tbcb_mech/atomistic
PY="$M/mm/envs/amb/bin/python -W ignore"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
O=interim/$1; mkdir -p $O
$PY scripts/a09_fep_analysis.py runs/fep --out $O/fep_interim.json > $O/a09.log 2>&1; echo "a09 exit $?" >> $O/a09.log
$PY scripts/a27_fep_convergence.py runs/fep --out $O/fep_convergence.json > $O/a27.log 2>&1; echo "a27 exit $?" >> $O/a27.log
ls -la $O
