#!/bin/bash
#SBATCH --job-name=CUTRUN_02_map
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --array=1-30%8
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Title:        Step 02 - alignment to GRCh38
# Input:        fastq_trimmed/<sample>_R{1,2}_trimmed.fastq.gz
# Output:       bam_files/<sample>_sorted.bam (+ .bai)
#               bam_files/qc/<sample>.flagstat.txt
#               logs/<sample>.bowtie2.log
# Depends on:   Bowtie2, SAMtools 1.11
# Notes:        Same Bowtie2 parameters as the CUT&Tag data:
#               --local --very-sensitive --no-unal --no-mixed --no-discordant
#               --phred33 -I 10 -X 700.
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"

if type module >/dev/null 2>&1; then
    module load bowtie2 2>/dev/null || true
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
fi
for tool in bowtie2 samtools; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found" >&2; exit 1; }
done
mkdir -p "${BAM_DIR}/qc"

PPN="${SLURM_CPUS_PER_TASK:-${THREADS}}"
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}" | tr -d '[:space:]')
[ -n "${SAMPLE}" ] || { echo "ERROR: no sample at index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

R1="${TRIM_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
R2="${TRIM_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
[ -f "${R1}" ] && [ -f "${R2}" ] || { echo "ERROR: trimmed FASTQ not found for ${SAMPLE}; run step 01" >&2; exit 1; }

OUT_BAM="${BAM_DIR}/${SAMPLE}_sorted.bam"
if [ -f "${OUT_BAM}" ]; then echo "Skipping ${SAMPLE}: BAM exists"; exit 0; fi

echo "[$(date)] Mapping ${SAMPLE}"
bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p "${PPN}" -x "${BOWTIE2_INDEX}" -1 "${R1}" -2 "${R2}" \
    2> "${LOG_DIR}/${SAMPLE}.bowtie2.log" \
  | samtools sort -@ "${PPN}" -O bam -o "${OUT_BAM}"

samtools index -@ "${PPN}" "${OUT_BAM}"
samtools flagstat "${OUT_BAM}" > "${BAM_DIR}/qc/${SAMPLE}.flagstat.txt"
echo "Done: ${SAMPLE} $(date)"
