#!/bin/bash
#SBATCH --job-name=cellbender
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-6
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 2b. Nucleus calling and ambient RNA removal, HuH6 (CellBender 0.4.0)
# FPR 0.01, 150 epochs. Expected nuclei / total droplets set per sample.

set -euo pipefail

CELLRANGER_DIR="cellranger"
OUT_DIR="cellbender"

# sample : expected nuclei : total droplets included
SAMPLES=(
    "HuH6_C70_1:3500:25000"
    "HuH6_C70_2:2500:20000"
    "HuH6_DMSO_1:3000:25000"
    "HuH6_DMSO_2:3500:25000"
    "HuH6_T5_1:2500:20000"
    "HuH6_T5_2:4500:25000"
)
IFS=':' read -r SAMPLE EXPECTED TOTAL <<< "${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}"

mkdir -p "${OUT_DIR}/${SAMPLE}"
cellbender remove-background \
    --input "${CELLRANGER_DIR}/${SAMPLE}/outs/raw_feature_bc_matrix.h5" \
    --output "${OUT_DIR}/${SAMPLE}/${SAMPLE}_cellbender.h5" \
    --expected-cells "${EXPECTED}" \
    --total-droplets-included "${TOTAL}" \
    --fpr 0.01 \
    --epochs 150 \
    --cpu-threads "${SLURM_CPUS_PER_TASK}"
