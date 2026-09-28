#!/bin/bash
#SBATCH --job-name=CUTRUN_07_windows
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# ============================================================
# Title:        Step 07 - summit windows per target
# Input:        results/peaks/<Condition>/<Target>/macs_q0.1_summits.bed
# Output:       results/windows/win_<target>_<width>.bed
#               results/windows/window_summary.tsv
# Depends on:   bedtools, SAMtools, Python 3
# Notes:        Summits from the three conditions are pooled per target.
#               Fixed windows are centred on summits (+/-250 bp for ATF2, Fos
#               and Jun; +/-1 kb for H3K4me3 and KDM5A), taking the strongest
#               summit first so that windows do not overlap. Windows that
#               overlap the blacklist are removed.
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"
export LC_ALL=C

if type module >/dev/null 2>&1; then
    module load bedtools 2>/dev/null || true
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
fi
for tool in bedtools samtools python3; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found" >&2; exit 1; }
done

TF_FLANK=250
BROAD_FLANK=1000
PEAKS="${RESULTS_DIR}/peaks"
OUT="${RESULTS_DIR}/windows"
TMP="${OUT}/.tmp"
rm -rf "${TMP}"; mkdir -p "${OUT}" "${TMP}"
trap 'rm -rf "${TMP}"' EXIT

REF_BAM=$(ls "${DEDUP_DIR}"/*_dedup.bam 2>/dev/null | head -1)
BAM_STYLE=$(samtools view -H "${REF_BAM}" | awk '/^@SQ/{sub(/^SN:/,"",$2); print $2; exit}')

# Main chromosomes only, in BAM naming.
main_chroms () {
    awk -v style="${BAM_STYLE}" 'BEGIN{OFS="\t"; want=(style ~ /^chr/)}
        NF>=2 && $1 !~ /^#/ { c=$1; sub(/^chr/,"",c)
            if(!(c ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/)) next
            $1=(want ? "chr"c : c); print }' "$1"
}
main_chroms "$(ensure_chrom_sizes)" | cut -f1,2 > "${TMP}/genome.sizes"
main_chroms "${BLACKLIST_BED3}" | cut -f1-3 | sort -k1,1 -k2,2n > "${TMP}/blacklist.bed"

# Pool summits (chrom, pos, pos+1, ., score) across conditions.
pool_summits () {
    local t="$1" c f
    : > "${TMP}/${t}.summits"
    for c in "${ALL_CONDITIONS[@]}"; do
        f="${PEAKS}/${c}/${t}/macs_q${MACS_QVALUE}_summits.bed"
        [ -s "${f}" ] || { echo "ERROR: missing ${f}; run step 06" >&2; exit 1; }
        main_chroms "${f}" | awk 'BEGIN{OFS="\t"} {print $1,$2,$2+1,".",$5+0}' >> "${TMP}/${t}.summits"
    done
}

# Place non-overlapping windows, strongest summit first; drop windows past chromosome ends.
place_windows () {
    local t="$1" flank="$2" out="$3"
    python3 - "${TMP}/${t}.summits" "${TMP}/genome.sizes" "${flank}" "${TMP}/${t}.placed" <<'PY'
import bisect, sys
from collections import defaultdict
summits, sizes_path, flank, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
sizes = {l.split()[0]: int(l.split()[1]) for l in open(sizes_path) if l.strip()}
best = {}
for l in open(summits):
    p = l.split()
    key, score = (p[0], int(p[1])), float(p[4])
    if key not in best or score > best[key]:
        best[key] = score
want = 2 * flank
taken, placed = defaultdict(list), []
for (chrom, pos), score in sorted(best.items(), key=lambda x: (-x[1], x[0][0], x[0][1])):
    t = taken[chrom]; i = bisect.bisect_left(t, pos)
    if (i > 0 and pos - t[i - 1] < want) or (i < len(t) and t[i] - pos < want):
        continue
    lo, hi = pos - flank, pos + flank
    if chrom not in sizes or lo < 0 or hi > sizes[chrom]:
        continue
    bisect.insort(t, pos); placed.append((chrom, lo, hi))
with open(out, "w") as fh:
    for chrom, lo, hi in sorted(placed):
        fh.write(f"{chrom}\t{lo}\t{hi}\n")
PY
    sort -k1,1 -k2,2n "${TMP}/${t}.placed" \
      | bedtools intersect -v -a - -b "${TMP}/blacklist.bed" \
      | awk -v p="${t}" 'BEGIN{OFS="\t"} {print $1,$2,$3,p"_"NR,".","+"}' > "${out}"
}

printf "window_set\tn_windows\twidth_bp\n" > "${OUT}/window_summary.tsv"
for T in "${MARKERS[@]}"; do
    lc=$(echo "${T}" | tr '[:upper:]' '[:lower:]')
    case "${T}" in H3K4me3|KDM5A) flank=${BROAD_FLANK}; tag="2kb" ;; *) flank=${TF_FLANK}; tag="500bp" ;; esac
    pool_summits "${T}"
    place_windows "${T}" "${flank}" "${OUT}/win_${lc}_${tag}.bed"
    printf "win_%s_%s\t%s\t%s\n" "${lc}" "${tag}" "$(wc -l < "${OUT}/win_${lc}_${tag}.bed" | tr -d ' ')" "$((2 * flank))" >> "${OUT}/window_summary.tsv"
done
cat "${OUT}/window_summary.tsv"
echo "Done: $(date)"
