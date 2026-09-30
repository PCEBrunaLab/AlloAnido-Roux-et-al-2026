#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-12
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# 1. Cell Ranger count
# snRNA-seq of SK-N-SH and HuH6, DMSO / KDM5-C70 / T-5224, 2 replicates each (12 libraries).
# Cell Ranger 8.0.0, reference refdata-gex-GRCh38-2020-A, Single Cell 3' v4.

set -euo pipefail

FASTQ_DIR="/path/to/fastq"
TRANSCRIPTOME="/path/to/refdata-gex-GRCh38-2020-A"
OUTPUT_DIR="cellranger"
SAMPLE_SHEET="sample_sheet.tsv"

# Row N of the sample sheet (after the header) for array task N
read -r SAMPLE_ID CELL_LINE TREATMENT REPLICATE CONDITION CONDITION_REP RUN LANE < \
    <(awk -v n="$((SLURM_ARRAY_TASK_ID + 1))" 'NR == n' "${SAMPLE_SHEET}")

mkdir -p "${OUTPUT_DIR}" && cd "${OUTPUT_DIR}"
echo "Cell Ranger: ${SAMPLE_ID} -> ${CONDITION_REP}"

cellranger count \
    --id="${CONDITION_REP}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${SAMPLE_ID}" \
    --create-bam=true \
    --localcores="${SLURM_CPUS_PER_TASK}" \
    --localmem=120
