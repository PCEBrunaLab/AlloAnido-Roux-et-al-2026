#!/bin/bash
#SBATCH --job-name=CUTRUN_09_promoters
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# ============================================================
# Title:        Step 09 - promoter intervals from SK-N-SH CUT&Tag H3K4me3
# Input:        ${CUTTAG_DIR}/BAM_dedup/SKNSH_H3K4Me3_Rep{1,2}_dedup.bam
# Output:       results/promoters/promoter.bed           (chr start end id score . summit)
#               results/promoters/promoter_summits.bed
#               results/promoters/promoter_summary.tsv
# Depends on:   MACS3 3.0.4, SAMtools, bedtools
# Notes:        The two CUT&Tag H3K4me3 replicates (processed with the same
#               trimming, alignment and deduplication) are merged and called
#               with the same MACS3 settings as the CUT&RUN data. Peaks on main
#               chromosomes that do not overlap the blacklist are the promoter
#               intervals. Chromosome names follow the CUT&RUN BAMs.
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"
export LC_ALL=C

MACS_BIN="$(find_macs)" || exit 1
if type module >/dev/null 2>&1; then
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
    module load bedtools 2>/dev/null || true
fi
for tool in samtools bedtools; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found" >&2; exit 1; }
done

P="${SLURM_CPUS_PER_TASK:-4}"
OUT="${RESULTS_DIR}/promoters"
TMP="${OUT}/.tmp"
mkdir -p "${OUT}" "${TMP}"
trap 'rm -rf "${TMP}"' EXIT

# Chromosome naming of the CUT&RUN BAMs.
REF_BAM=$(ls "${MERGED_BAM_DIR}"/*_merged_sorted.bam 2>/dev/null | head -1)
[ -n "${REF_BAM}" ] || { echo "ERROR: no merged CUT&RUN BAM; run step 05" >&2; exit 1; }
STYLE=$(samtools view -H "${REF_BAM}" | awk '/^@SQ/{sub(/^SN:/,"",$2); print $2; exit}')
to_cr_style () {
    awk -v want="$([[ "${STYLE}" == chr* ]] && echo 1 || echo 0)" 'BEGIN{OFS="\t"} NF>=3 && $1 !~ /^#/ {
        c=$1; sub(/^chr/,"",c)
        if(!(c ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/)) next
        $1=(want==1 ? "chr"c : c); print }' "$1"
}
to_cr_style "${BLACKLIST_BED3}" | cut -f1-3 | sort -k1,1 -k2,2n > "${TMP}/blacklist.bed"

# Merge the two CUT&Tag replicates.
REPS=()
for R in 1 2; do
    b="${CUTTAG_DIR}/BAM_dedup/SKNSH_H3K4Me3_Rep${R}_dedup.bam"
    [ -s "${b}" ] || { echo "ERROR: missing ${b}" >&2; exit 1; }
    REPS+=("${b}")
done
MERGED="${OUT}/SKNSH_H3K4me3_merged_sorted.bam"
if [ ! -s "${MERGED}" ]; then
    samtools merge -@ "${P}" -f "${TMP}/merged.bam" "${REPS[@]}"
    samtools sort -@ "${P}" -o "${MERGED}" "${TMP}/merged.bam"
    samtools index "${MERGED}"
fi

# Peak calling, same settings as step 06.
NAME="macs_q${MACS_QVALUE}"
"${MACS_BIN}" callpeak -t "${MERGED}" -f "${PEAK_FORMAT}" -g "${GENOME_CODE}" -q "${MACS_QVALUE}" \
    --nomodel --keep-dup all -n "${NAME}" --outdir "${TMP}/" 2> "${OUT}/${NAME}.log"
[ -s "${TMP}/${NAME}_peaks.narrowPeak" ] || { echo "ERROR: no peaks written" >&2; exit 1; }

# narrowPeak -> chr start end id score(-log10 q) . absolute_summit
to_cr_style "${TMP}/${NAME}_peaks.narrowPeak" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,".",$9,".",$2+$10}' \
  | sort -k1,1 -k2,2n \
  | bedtools intersect -v -a - -b "${TMP}/blacklist.bed" \
  | awk 'BEGIN{OFS="\t"} {$4="promoter_"NR; print}' > "${OUT}/promoter.bed"
awk 'BEGIN{OFS="\t"} {print $1,$7,$7+1,$4,$5,$6}' "${OUT}/promoter.bed" > "${OUT}/promoter_summits.bed"

awk '{w=$3-$2; s+=w; print w}' "${OUT}/promoter.bed" | sort -n \
  | awk '{w[NR]=$1; s+=$1} END{n=NR; m=(n%2 ? w[(n+1)/2] : (w[n/2]+w[n/2+1])/2);
         printf "n\tmean_bp\tmedian_bp\n%d\t%.0f\t%.0f\n", n, s/n, m}' > "${OUT}/promoter_summary.tsv"
cat "${OUT}/promoter_summary.tsv"
echo "Done: $(date)"
