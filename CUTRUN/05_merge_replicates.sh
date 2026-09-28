#!/bin/bash
#SBATCH --job-name=CUTRUN_05_merge
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# ============================================================
# Title:        Step 05 - pool replicates per condition and target
# Input:        bam_files_dedup/<sample>_dedup.bam
# Output:       merged_bam_files/<Condition>_<Target>_merged_sorted.bam (+ .bai)
#               merged_bam_files/<Condition>_<Target>_merged.flagstat.txt
# Depends on:   SAMtools 1.11
# Notes:        The two deduplicated replicates are merged
#               (3 conditions x 5 targets = 15 BAMs).
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"

if type module >/dev/null 2>&1; then
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
fi
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }

N=0
for condition in "${ALL_CONDITIONS[@]}"; do
    for target in "${MARKERS[@]}"; do
        GROUP="${condition}_${target}"
        OUT="${MERGED_BAM_DIR}/${GROUP}_merged_sorted.bam"
        IN=()
        for rep in "${REPLICATES[@]}"; do
            sid="${sample_map[${condition}_${target}_${rep}]:-}"
            bam="${DEDUP_DIR}/${sid}_dedup.bam"
            [ -n "${sid}" ] && [ -f "${bam}" ] || { echo "ERROR: missing ${GROUP} replicate ${rep}" >&2; exit 1; }
            IN+=("${bam}")
        done
        if [ ! -f "${OUT}" ]; then
            echo "Merging ${GROUP}"
            samtools merge -@ "${THREADS}" -f "${OUT}.unsorted.bam" "${IN[@]}"
            samtools sort -@ "${THREADS}" -O bam -o "${OUT}" "${OUT}.unsorted.bam"
            rm -f "${OUT}.unsorted.bam"
            samtools index -@ "${THREADS}" "${OUT}"
            samtools flagstat "${OUT}" > "${MERGED_BAM_DIR}/${GROUP}_merged.flagstat.txt"
        fi
        N=$((N + 1))
    done
done

echo "Merged BAMs: ${N}"
echo "Done: $(date)"
