#!/bin/bash
#SBATCH --job-name=CUTRUN_08_tracks
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# ============================================================
# Title:        Step 08 - CPM coverage tracks and occupancy heatmaps
# Input:        merged_bam_files/<Condition>_<Target>_merged_sorted.bam
#               results/windows/win_<target>_<width>.bed
# Output:       results/bigwig/<Condition>_<Target>.cpm.bw
#               results/heatmaps/<Target>_profile_heatmap.pdf
# Depends on:   deepTools 3.5.5, SAMtools
# Notes:        bamCoverage --normalizeUsing CPM --binSize 50 --extendReads,
#               blacklist excluded. One computeMatrix per target over all
#               conditions (+/-3 kb, 50 bp bins), so rows are identical
#               across conditions.
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"

activate_deeptools computeMatrix || exit 1
if type module >/dev/null 2>&1; then
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
fi
for tool in bamCoverage computeMatrix plotHeatmap samtools; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found" >&2; exit 1; }
done

P="${SLURM_CPUS_PER_TASK:-${THREADS}}"
BWDIR="${RESULTS_DIR}/bigwig"
HMDIR="${RESULTS_DIR}/heatmaps"
mkdir -p "${BWDIR}" "${HMDIR}/mat"

for T in "${MARKERS[@]}"; do
    lc=$(echo "${T}" | tr '[:upper:]' '[:lower:]')
    case "${T}" in H3K4me3|KDM5A) bed="${RESULTS_DIR}/windows/win_${lc}_2kb.bed" ;; *) bed="${RESULTS_DIR}/windows/win_${lc}_500bp.bed" ;; esac
    [ -s "${bed}" ] || { echo "ERROR: missing ${bed}; run step 07" >&2; exit 1; }
    NREG=$(wc -l < "${bed}" | tr -d ' ')
    echo "==== ${T} (${NREG} windows) ===="

    BWS=()
    for C in "${ALL_CONDITIONS[@]}"; do
        bam="${MERGED_BAM_DIR}/${C}_${T}_merged_sorted.bam"
        bw="${BWDIR}/${C}_${T}.cpm.bw"
        [ -s "${bam}" ] || { echo "ERROR: missing ${bam}; run step 05" >&2; exit 1; }
        [ -s "${bw}" ] || bamCoverage -b "${bam}" -o "${bw}" \
            --normalizeUsing CPM --binSize "${HEATMAP_BINSIZE}" --extendReads \
            --blackListFileName "${BLACKLIST_BED3}" -p "${P}"
        BWS+=("${bw}")
    done

    mat="${HMDIR}/mat/${T}.mat.gz"
    computeMatrix reference-point --referencePoint center -R "${bed}" -S "${BWS[@]}" \
        -b "${HEATMAP_DISTANCE}" -a "${HEATMAP_DISTANCE}" --binSize "${HEATMAP_BINSIZE}" \
        --missingDataAsZero --samplesLabel "${CONDITION_LABELS[@]}" -o "${mat}" -p "${P}"

    ZMAX=(--zMin 0 --zMax auto)
    [ "${T}" = "H3K4me3" ] && ZMAX=(--zMin 0 --zMax 3)
    plotHeatmap -m "${mat}" -out "${HMDIR}/${T}_profile_heatmap.pdf" \
        --colorMap RdYlBu_r --whatToShow 'plot, heatmap and colorbar' \
        --heatmapHeight 14 --heatmapWidth 4 --sortUsing mean --sortRegions descend \
        --refPointLabel summit --plotTitle "${T} (CPM, n = ${NREG})" \
        --xAxisLabel "bp" --yAxisLabel "" --regionsLabel "" --legendLocation none \
        "${ZMAX[@]}" --dpi 200
done
echo "Done: $(date)"
