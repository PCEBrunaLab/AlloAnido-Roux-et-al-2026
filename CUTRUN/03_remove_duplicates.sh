#!/bin/bash
#SBATCH --job-name=CUTRUN_03_dedup
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# ============================================================
# Title:        Step 03 - duplicate removal
# Input:        bam_files/<sample>_sorted.bam
# Output:       bam_files_dedup/<sample>_dedup.bam (+ .bai)
#               bam_files_dedup/metrics/<sample>.dedup_metrics.txt
#               qc/qc_summary_postdedup.txt
# Depends on:   Picard 2.23.8, SAMtools 1.11
# Notes:        MarkDuplicates with REMOVE_DUPLICATES=true.
#               No mapping-quality filter is applied.
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"

if type module >/dev/null 2>&1; then
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
    module load picard 2>/dev/null || module load picard-tools 2>/dev/null || true
fi
for tool in samtools picard; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found" >&2; exit 1; }
done

METRICS_DIR="${DEDUP_DIR}/metrics"
mkdir -p "${METRICS_DIR}"
SUMMARY="${DEDUP_QC_DIR}/qc_summary_postdedup.txt"
printf "Sample\tTotal_Reads\tMapped\tDuplicates_Removed\tPercent_Duplicates\tUnique_Fragments\tCondition\tTarget\tReplicate\n" > "${SUMMARY}"

while IFS= read -r SAMPLE || [ -n "${SAMPLE}" ]; do
    SAMPLE=$(echo "${SAMPLE}" | tr -d '[:space:]\r')
    [ -z "${SAMPLE}" ] && continue
    IN="${BAM_DIR}/${SAMPLE}_sorted.bam"
    OUT="${DEDUP_DIR}/${SAMPLE}_dedup.bam"
    MET="${METRICS_DIR}/${SAMPLE}.dedup_metrics.txt"
    [ -f "${IN}" ] || { echo "Missing ${IN}, skipping"; continue; }

    if [ ! -f "${OUT}" ]; then
        echo "[$(date)] Deduplicating ${SAMPLE}"
        picard MarkDuplicates I="${IN}" O="${OUT}" M="${MET}" \
            REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp
        samtools index "${OUT}"
    fi

    TOTAL=$(samtools view -c "${OUT}")
    MAPPED=$(samtools view -c -F 4 "${OUT}")
    UNIQUE=$(samtools view -c -f 2 "${OUT}" | awk '{print int($1/2)}')
    DUPS=$(grep -A 1 "LIBRARY" "${MET}" | tail -1 | cut -f 8)
    PCT=$(grep -A 1 "LIBRARY" "${MET}" | tail -1 | cut -f 9 | awk '{printf "%.1f", $1*100}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${SAMPLE}" "${TOTAL}" "${MAPPED}" "${DUPS}" "${PCT}" "${UNIQUE}" \
        "${sample_condition[$SAMPLE]:-NA}" "${sample_target[$SAMPLE]:-NA}" "${sample_replicate[$SAMPLE]:-NA}" >> "${SUMMARY}"
done < "${SAMPLE_LIST}"

echo "Done: $(date)"
