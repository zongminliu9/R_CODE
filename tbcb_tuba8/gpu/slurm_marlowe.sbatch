#!/bin/bash
#SBATCH --job-name=tbcb-fep
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# Marlowe's H100s roughly halve the wall times quoted in README.md, so the 48 h
# Sherlock window is cut to 24 h here. Set --partition / --account to your allocation;
# `sinfo -s` lists what you can use.
set -euo pipefail
JOB="${1:?usage: sbatch slurm_marlowe.sbatch <A1|A2|A3|B1|B2|C1> [extra args]}"; shift || true
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate tbcb
mkdir -p logs out
export OPENMM_DEFAULT_PLATFORM=CUDA
# H100 benefits from TF32 off for FEP reproducibility
export OPENMM_CUDA_PRECISION=mixed
exec bash slurm_sherlock.sbatch "$JOB" "$@"
