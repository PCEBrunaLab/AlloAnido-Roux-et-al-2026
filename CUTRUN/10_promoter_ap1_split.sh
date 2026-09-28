#!/bin/bash
#SBATCH --job-name=CUTRUN_10_ap1_split
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=48G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# ============================================================
# Title:        Step 10 - AP-1-high and AP-1-low promoters within the top
#               H3K4me3 quartile
# Input:        results/promoters/promoter.bed
#               results/bigwig/<Condition>_<Target>.cpm.bw
# Output:       results/promoter_ap1_split/results.tsv
#               results/promoter_ap1_split/<Factor>/summary.tsv
#               results/promoter_ap1_split/<Factor>/regions/{high,low}.bed
#               results/promoter_ap1_split/<Factor>/{KDM5A,H3K4me3,<Factor>}_heatmap.pdf
# Depends on:   deepTools 3.5.5, Python 3 (numpy)
# Notes:        1. Mean CPM of H3K4me3, KDM5A, ATF2, Fos and Jun over each
#                  promoter interval (multiBigwigSummary).
#               2. Ranking scores are averaged over SPLIT_CONDS.
#               3. The top 25% of promoters by H3K4me3 are kept.
#               4. Within them, the top and bottom 25% by ATF2, Fos or Jun
#                  define AP-1-high and AP-1-low.
#               5. KDM5A and H3K4me3 are compared per condition as the ratio
#                  of median CPM (high / low).
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
[ -f "${SCRIPT_DIR}/00_config.sh" ] || SCRIPT_DIR="${SLURM_SUBMIT_DIR:-.}"
source "${SCRIPT_DIR}/00_config.sh"
export LC_ALL=C

activate_deeptools computeMatrix || exit 1
for tool in computeMatrix plotHeatmap multiBigwigSummary python3; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not found" >&2; exit 1; }
done

P="${SLURM_CPUS_PER_TASK:-${THREADS}}"
H3_PCT=25
AP1_PCT=25
SPLIT_CONDS="${SPLIT_CONDS:-${ALL_CONDITIONS[*]}}"
ALL_CONDS="${ALL_CONDITIONS[*]}"
FACTORS="ATF2 Fos Jun"
ZQ=0.98                                   # colour scale: 98th percentile per panel

PROM="${RESULTS_DIR}/promoters/promoter.bed"
BW="${RESULTS_DIR}/bigwig"
OUT="${RESULTS_DIR}/promoter_ap1_split"
TMP="${OUT}/.tmp"
[ -s "${PROM}" ] || { echo "ERROR: ${PROM} missing; run step 09" >&2; exit 1; }
for C in ${ALL_CONDS}; do for T in "${MARKERS[@]}"; do
    [ -s "${BW}/${C}_${T}.cpm.bw" ] || { echo "ERROR: missing ${BW}/${C}_${T}.cpm.bw; run step 08" >&2; exit 1; }
done; done
rm -rf "${TMP}"; mkdir -p "${OUT}" "${TMP}"
trap 'rm -rf "${TMP}"' EXIT
echo "promoters -> top ${H3_PCT}% H3K4me3 -> AP-1 top/bottom ${AP1_PCT}%   ranking conditions: ${SPLIT_CONDS}"

# ---- 1. mean CPM per promoter ----------------------------------------------
cut -f1-3 "${PROM}" | sort -k1,1 -k2,2n > "${TMP}/promoters.bed3"
awk 'BEGIN{OFS="\t"} {print $1":"$2":"$3, $7}' "${PROM}" | sort -k1,1 > "${TMP}/summit.txt"

BWS=(); LABS=()
for C in ${SPLIT_CONDS}; do for T in ATF2 Fos Jun; do BWS+=("${BW}/${C}_${T}.cpm.bw"); LABS+=("${T}_${C}"); done; done
for C in ${ALL_CONDS}; do
    BWS+=("${BW}/${C}_H3K4me3.cpm.bw" "${BW}/${C}_KDM5A.cpm.bw"); LABS+=("H3K4me3_${C}" "KDM5A_${C}")
done
multiBigwigSummary BED-file -b "${BWS[@]}" --labels "${LABS[@]}" --BED "${TMP}/promoters.bed3" \
    -p "${P}" -o "${TMP}/signal.npz" --outRawCounts "${TMP}/signal.tab" >/dev/null 2>&1
[ -s "${TMP}/signal.tab" ] || { echo "ERROR: multiBigwigSummary produced no output" >&2; exit 1; }

# ---- 2. selection and comparison -------------------------------------------
python3 - "${TMP}/signal.tab" "${TMP}/summit.txt" "${OUT}" "${H3_PCT}" "${AP1_PCT}" \
    "${SPLIT_CONDS}" "${FACTORS}" "${ALL_CONDS}" <<'PY'
import sys, os, math, statistics as st
raw, summit_p, out = sys.argv[1:4]
h3_pct, ap1_pct = float(sys.argv[4]) / 100, float(sys.argv[5]) / 100
SPLIT, FACTORS, ALL = sys.argv[6].split(), sys.argv[7].split(), sys.argv[8].split()

rows, hdr = [], None
for ln in open(raw):
    if hdr is None:
        hdr = [c.strip().strip("'\"") for c in ln.lstrip("#").rstrip("\n").split("\t")]
        hdr[0] = "chr"
        continue
    p = ln.rstrip("\n").split("\t")
    if len(p) >= len(hdr):
        rows.append(dict(zip(hdr, p)))

def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")
def mean(v):
    v = [x for x in v if not math.isnan(x)]
    return sum(v) / len(v) if v else float("nan")
def med(v):
    v = [x for x in v if not math.isnan(x)]
    return st.median(v) if v else float("nan")
def key(d):
    return f"{d['chr']}:{d['start']}:{d['end']}"

summit = dict(l.split() for l in open(summit_p))

pool = []
for d in rows:
    h = mean([num(d.get(f"H3K4me3_{C}")) for C in SPLIT])
    if not math.isnan(h):
        d["_h3"] = h
        pool.append(d)
pool.sort(key=lambda d: -d["_h3"])
top_h3 = pool[:max(1, int(round(len(pool) * h3_pct)))]
print(f"promoters: {len(rows)}   top {int(h3_pct*100)}% H3K4me3: {len(top_h3)}")

res = open(os.path.join(out, "results.tsv"), "w")
res.write("factor\tn_promoters\tn_top_H3K4me3\tn_high\tn_low\t" +
          "\t".join(f"KDM5A_ratio_{C}" for C in ALL) + "\tH3K4me3_ratio\n")
for F in FACTORS:
    scored = []
    for d in top_h3:
        a = mean([num(d.get(f"{F}_{C}")) for C in SPLIT])
        if not math.isnan(a):
            d = dict(d); d["_ap1"] = a; scored.append(d)
    scored.sort(key=lambda d: (-d["_ap1"], key(d)))
    k = max(1, int(len(scored) * ap1_pct))
    hi, lo = scored[:k], scored[-k:]
    def kd(d):                                                           # heatmap row order: mean KDM5A, descending
        m = mean([num(d.get(f"KDM5A_{C}")) for C in ALL])
        return 0.0 if math.isnan(m) else -m
    hi.sort(key=kd); lo.sort(key=kd)

    od = os.path.join(out, F)
    os.makedirs(os.path.join(od, "regions"), exist_ok=True)
    for tag, grp in (("high", hi), ("low", lo)):
        with open(os.path.join(od, "regions", f"{tag}.bed"), "w") as fh:
            for d in grp:
                s = int(float(summit.get(key(d), (int(d["start"]) + int(d["end"])) // 2)))
                fh.write(f"{d['chr']}\t{s}\t{s + 1}\n")

    ratios = []
    with open(os.path.join(od, "summary.tsv"), "w") as fh:
        fh.write("condition\tn_high\tn_low\tKDM5A_high\tKDM5A_low\tKDM5A_ratio\tH3K4me3_high\tH3K4me3_low\tH3K4me3_ratio\n")
        for C in ALL:
            kh, kl = med([num(d.get(f"KDM5A_{C}")) for d in hi]), med([num(d.get(f"KDM5A_{C}")) for d in lo])
            hh, hl = med([num(d.get(f"H3K4me3_{C}")) for d in hi]), med([num(d.get(f"H3K4me3_{C}")) for d in lo])
            kr, hr = kh / kl if kl else float("nan"), hh / hl if hl else float("nan")
            ratios.append(kr)
            fh.write(f"{C}\t{len(hi)}\t{len(lo)}\t{kh:.4g}\t{kl:.4g}\t{kr:.3f}\t{hh:.4g}\t{hl:.4g}\t{hr:.3f}\n")
    h3r = med([d["_h3"] for d in hi]) / med([d["_h3"] for d in lo])
    res.write(f"{F}\t{len(rows)}\t{len(top_h3)}\t{len(hi)}\t{len(lo)}\t" +
              "\t".join(f"{r:.3f}" for r in ratios) + f"\t{h3r:.3f}\n")
    print(f"{F}: {len(hi)} vs {len(lo)}   KDM5A ratio " +
          "  ".join(f"{C} {r:.2f}" for C, r in zip(ALL, ratios)) + f"   H3K4me3 ratio {h3r:.2f}")
res.close()
PY

# ---- 3. heatmaps -------------------------------------------------------------
zcap () {
    python3 - "$1" "${ZQ}" <<'PY'
import sys, gzip, numpy as np
v = []
with gzip.open(sys.argv[1], "rt") as fh:
    for l in fh:
        if not l.startswith("@"):
            v.append(np.array(l.rstrip("\n").split("\t")[6:], dtype=float))
v = np.concatenate(v); v = v[np.isfinite(v)]
print("%.4g" % (np.quantile(v, float(sys.argv[2])) if v.size else 1.0))
PY
}
heatmap () {                              # $1 dir  $2 mark  $3 title
    local od="$1" mark="$2" title="$3" bws=() mat
    for C in ${ALL_CONDS}; do bws+=("${BW}/${C}_${mark}.cpm.bw"); done
    mat="${TMP}/$(basename "${od}")_${mark}.mat.gz"
    computeMatrix reference-point --referencePoint center -b "${HEATMAP_DISTANCE}" -a "${HEATMAP_DISTANCE}" \
        --binSize "${HEATMAP_BINSIZE}" -R "${od}/regions/high.bed" "${od}/regions/low.bed" -S "${bws[@]}" \
        --samplesLabel ${ALL_CONDS} --sortRegions keep --missingDataAsZero -p "${P}" -o "${mat}" >/dev/null 2>&1
    plotHeatmap -m "${mat}" -o "${od}/${mark}_heatmap.pdf" --sortRegions keep \
        --colorMap RdYlBu_r --zMin 0 --zMax "$(zcap "${mat}")" --yMin 0 \
        --heatmapHeight 12 --heatmapWidth 4 --refPointLabel "summit" --xAxisLabel "" --yAxisLabel "CPM" \
        --regionsLabel "${LH}" "${LL}" --legendLocation none --plotTitle "${title}" \
        --plotFileFormat pdf --dpi 200 >/dev/null 2>&1
    echo "  -> ${od}/${mark}_heatmap.pdf"
}
for F in ${FACTORS}; do
    od="${OUT}/${F}"
    LH="${F}-high (n=$(wc -l < "${od}/regions/high.bed" | tr -d ' '))"
    LL="${F}-low (n=$(wc -l < "${od}/regions/low.bed" | tr -d ' '))"
    heatmap "${od}" KDM5A   "KDM5A at top-${H3_PCT}% H3K4me3 promoters, split by ${F}"
    heatmap "${od}" H3K4me3 "H3K4me3 at top-${H3_PCT}% H3K4me3 promoters, split by ${F}"
    heatmap "${od}" "${F}"  "${F} at top-${H3_PCT}% H3K4me3 promoters, split by ${F}"
done

column -t -s $'\t' "${OUT}/results.tsv"
echo "Done: $(date)"
