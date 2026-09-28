#!/bin/bash
#SBATCH --job-name=CUTRUN_01_trim
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-30%8
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Title:        Step 01 - adapter and poly-G trimming
# Input:        raw_fastq/<sample>_R{1,2}_001.fastq.gz
# Output:       fastq_trimmed/<sample>_R{1,2}_trimmed.fastq.gz
#               fastq_trimmed/reports/<sample>.cutadapt.txt
# Depends on:   Cutadapt 4.9
# Notes:        TruSeq adapters; --nextseq-trim=20 -O 1 -q 0 -m 20.
#               One array task per line of sample_list.txt.
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"

if ! command -v cutadapt >/dev/null 2>&1 && type module >/dev/null 2>&1; then
    module load cutadapt 2>/dev/null || true
fi
command -v cutadapt >/dev/null 2>&1 || { echo "ERROR: cutadapt not found" >&2; exit 1; }

ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
REPORT_DIR="${TRIM_DIR}/reports"
mkdir -p "${REPORT_DIR}"

PPN="${SLURM_CPUS_PER_TASK:-${THREADS}}"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}" | tr -d '[:space:]')
[ -n "${SAMPLE}" ] || { echo "ERROR: no sample at index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

R1="${FASTQ_RUN_DIR}/${SAMPLE}_R1_001.fastq.gz"
R2="${FASTQ_RUN_DIR}/${SAMPLE}_R2_001.fastq.gz"
[ -f "${R1}" ] && [ -f "${R2}" ] || { echo "ERROR: FASTQ not found for ${SAMPLE}" >&2; exit 1; }

OUT_R1="${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
OUT_R2="${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
if [ -f "${OUT_R1}" ] && [ -f "${OUT_R2}" ]; then
    echo "Skipping ${SAMPLE}: trimmed FASTQ exists"; exit 0
fi

echo "[$(date)] Trimming ${SAMPLE}"
cutadapt \
    -a "${ADAPTER_R1}" -A "${ADAPTER_R2}" \
    -O 1 --nextseq-trim=20 -q 0 -m "${MIN_READ_LENGTH}" \
    -j "${PPN}" \
    -o "${OUT_R1}" -p "${OUT_R2}" \
    "${R1}" "${R2}" > "${REPORT_DIR}/${SAMPLE}.cutadapt.txt"

echo "Done: ${SAMPLE} $(date)"
