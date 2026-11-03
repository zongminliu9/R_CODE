#!/bin/bash
#SBATCH --job-name=tbcb-fep
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# Sherlock: add '-C GPU_MEM:40GB' or '-p owners' as your allocation allows.
# A 70k-atom ff19SB/OPC system fits comfortably in 16 GB, so any modern card works.
set -euo pipefail
JOB="${1:?usage: sbatch slurm_sherlock.sbatch <A1|A2|A3|B1|B2|C1> [extra args]}"; shift || true
module load cuda/12.4
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate tbcb
mkdir -p logs out
export OPENMM_DEFAULT_PLATFORM=CUDA
echo "host=$(hostname)  job=$JOB  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader)"
case "$JOB" in
  A1) python 10_fep_mutation.py --system systems/tuba1a_solv.xml --mutation I114S \
        --legs folded --replicas 3 --lambdas 20 --ns-per-window 10 --out out/A1 "$@" ;;
  A2) python 10_fep_mutation.py --system systems/fragment_I114.xml --mutation I114S \
        --legs reference --replicas 3 --lambdas 20 --ns-per-window 20 --out out/A2 "$@" ;;
  A3) python 10_fep_mutation.py --system systems/tuba1a_solv.xml \
        --mutation "${2:-I110V}" --legs folded --replicas 3 --lambdas 20 \
        --ns-per-window 10 --out "out/A3_${2:-I110V}" "$@" ;;
  B1) python 20_rest2_loop.py --pair tuba8 tuba8_loop1a --replicas 16 \
        --tmax 450 --ns 200 --out out/B1 "$@" ;;
  B2) python 20_rest2_loop.py --pair tuba1a tuba1a_loop8 --replicas 16 \
        --tmax 450 --ns 200 --out out/B2 "$@" ;;
  C1) python 10_fep_mutation.py --system systems/tbcb_tuba1a_solv.xml --mutation Y451F \
        --legs complex free --replicas 3 --lambdas 12 --ns-per-window 10 --out out/C1 "$@" ;;
  *) echo "unknown job $JOB" >&2; exit 2 ;;
esac
