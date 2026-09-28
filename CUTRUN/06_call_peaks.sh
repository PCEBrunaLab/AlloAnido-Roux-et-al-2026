#!/bin/bash
#SBATCH --job-name=CUTRUN_06_peaks
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --array=1-15
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

# ============================================================
# Title:        Step 06 - peak calling on pooled BAMs
# Input:        merged_bam_files/<Condition>_<Target>_merged_sorted.bam
# Output:       results/peaks/<Condition>/<Target>/macs_q0.1_peaks.narrowPeak.bed
#               results/peaks/<Condition>/<Target>/macs_q0.1_summits.bed
# Depends on:   MACS3 3.0.4, SAMtools
# Notes:        callpeak -f BAMPE -g hs -q 0.1 --nomodel --keep-dup all,
#               no control, local lambda. One array task per condition x target
#               (3 x 5 = 15).
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"

MACS_BIN="$(find_macs)" || exit 1
echo "Peak caller: ${MACS_BIN} ($(${MACS_BIN} --version 2>&1 | tr -d '\n'))"

N_COND=${#ALL_CONDITIONS[@]}
IDX="${SLURM_ARRAY_TASK_ID:-1}"
TARGET="${MARKERS[$(( (IDX - 1) / N_COND ))]}"
COND="${ALL_CONDITIONS[$(( (IDX - 1) % N_COND ))]}"

BAM="${MERGED_BAM_DIR}/${COND}_${TARGET}_merged_sorted.bam"
[ -s "${BAM}" ] || { echo "ERROR: missing ${BAM}; run step 05" >&2; exit 1; }

OUT="${RESULTS_DIR}/peaks/${COND}/${TARGET}"
NAME="macs_q${MACS_QVALUE}"
mkdir -p "${OUT}"

echo "=== ${COND} ${TARGET} ==="
"${MACS_BIN}" callpeak \
    -t "${BAM}" -f "${PEAK_FORMAT}" -g "${GENOME_CODE}" -q "${MACS_QVALUE}" \
    --nomodel --keep-dup all \
    -n "${NAME}" --outdir "${OUT}/" 2> "${OUT}/${NAME}.log"

[ -s "${OUT}/${NAME}_peaks.narrowPeak" ] || { echo "ERROR: no peaks written" >&2; tail -20 "${OUT}/${NAME}.log" >&2; exit 1; }
mv "${OUT}/${NAME}_peaks.narrowPeak" "${OUT}/${NAME}_peaks.narrowPeak.bed"
echo "${COND} ${TARGET}: $(wc -l < "${OUT}/${NAME}_peaks.narrowPeak.bed") peaks"
echo "Done: $(date)"
