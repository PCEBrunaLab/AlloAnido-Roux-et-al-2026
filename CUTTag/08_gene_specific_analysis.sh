#!/bin/bash
#SBATCH --job-name=CUTTag_GeneSpecific
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 08: Gene-Specific Peak Analysis
# Analyzes signal at specific gene sets (e.g., target genes)
# Creates TSS-centered and peak-centered heatmaps per gene list
#
# Input:  BigWig files (Step 06), gene list CSV, gene annotations
# Output: Gene-specific heatmaps, peak assignments, active gene sets
#==================================================================================

source $(dirname "$0")/00_config.sh

OUT_DIR="${HEATMAP_DIR}/gene_specific_peaks"
mkdir -p ${OUT_DIR}

HOMER_DIST=7500

echo "Gene-Specific Peak and TSS Analysis"
echo "Peak-to-gene assignment window: ${PEAK_GENE_DISTANCE} bp"

# Parse CSV gene list and create per-condition gene files
HEADER=$(head -n1 "${GENE_LIST_CSV}" | sed 's/^\xEF\xBB\xBF//' | tr -d '\r')
IFS=',' read -ra GENE_CONDITIONS <<< "${HEADER}"

declare -A GENE_LISTS

for gc in "${GENE_CONDITIONS[@]}"; do
    gc=$(echo "$gc" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//')
    [ -z "$gc" ] && continue

    # Find column index
    col_idx=1
    for i in "${!GENE_CONDITIONS[@]}"; do
        cond=$(echo "${GENE_CONDITIONS[i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//')
        if [ "$cond" = "$gc" ]; then
            col_idx=$((i+1))
            break
        fi
    done

    GENE_FILE="${OUT_DIR}/${gc}_genes.txt"
    tail -n +2 "${GENE_LIST_CSV}" | cut -d',' -f"${col_idx}" | \
        sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//' | \
        grep -v -E '^\s*$' | sort -u > "${GENE_FILE}"

    NUM=$(wc -l < "${GENE_FILE}" 2>/dev/null || echo "0")
    echo "  ${gc}: ${NUM} genes"
    [ "$NUM" -gt 0 ] && GENE_LISTS[${gc}]="${GENE_FILE}"
done

# Process each condition's gene list
for gc in "${!GENE_LISTS[@]}"; do
    GENE_FILE="${GENE_LISTS[${gc}]}"
    GC_DIR="${OUT_DIR}/${gc}"
    mkdir -p "${GC_DIR}"

    # Create gene subset BED and TSS BED
    SUBSET_BED="${GC_DIR}/genes_subset.bed"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${GENE_FILE}" "${GENES_BED}" | sort -k1,1 -k2,2n | uniq > "${SUBSET_BED}"

    NUM_GENES=$(wc -l < "${SUBSET_BED}" 2>/dev/null || echo "0")
    [ "$NUM_GENES" -eq 0 ] && continue

    TSS_BED="${GC_DIR}/tss.bed"
    awk 'BEGIN{OFS="\t"} {
        if(NF<6) next;
        if($6=="+") tss=$2; else if($6=="-") tss=$3-1; else next;
        if(tss<0) tss=0;
        print $1, tss, tss+1, $4"_TSS", 0, $6
    }' "${SUBSET_BED}" | sort -k1,1 -k2,2n | uniq > "${TSS_BED}"

    for marker in ${MARKERS[@]}; do
        # Build BigWig list
        BW_FILES=""
        LABELS=""
        for i in "${!ALL_CONDITIONS[@]}"; do
            bw="${BIGWIG_DIR}/${ORGANISM}/${ALL_CONDITIONS[$i]}/${marker}/normalized.bw"
            [ -f "$bw" ] && BW_FILES="${BW_FILES} ${bw}" && LABELS="${LABELS} ${CONDITION_LABELS[$i]}"
        done
        [ -z "$BW_FILES" ] && continue
        read -ra LABEL_ARRAY <<< "${LABELS}"

        # TSS heatmap
        computeMatrix reference-point \
            -S ${BW_FILES} -R "${TSS_BED}" \
            --referencePoint TSS -a 3000 -b 3000 \
            --skipZeros --missingDataAsZero \
            -o "${GC_DIR}/${marker}_TSS_matrix.gz" \
            --samplesLabel ${LABEL_ARRAY[@]} -p ${THREADS}

        plotHeatmap -m "${GC_DIR}/${marker}_TSS_matrix.gz" \
            -out "${GC_DIR}/${marker}_${gc}_TSS_heatmap.pdf" \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 15 --heatmapWidth 10 \
            --plotTitle "${marker} at ${gc} Gene TSS" \
            --refPointLabel "TSS" --sortUsing mean --sortRegions descend

        # Peak-centered analysis
        PEAK_FILE="${PROJ_DIR}/Homer/mergePeaks/${ORGANISM}/${TIME_CONTROL}_control/${marker}/${HOMER_DIST}/mergedPeakFile.bed"
        [ ! -f "$PEAK_FILE" ] && continue

        PEAK_CLEAN="${GC_DIR}/${marker}_peaks_clean.bed"
        awk 'BEGIN{OFS="\t"} {
            if(NF<3) next; if($2<0) $2=0; if($3<=$2) $3=$2+1;
            if($1 !~ /^chr/) $1="chr"$1;
            print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR)
        }' "${PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"

        FILTERED="${GC_DIR}/${marker}_gene_associated_peaks.bed"
        bedtools window -w ${PEAK_GENE_DISTANCE} -a "${PEAK_CLEAN}" -b "${SUBSET_BED}" | \
            awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_near_"$8}' | sort -u > "${FILTERED}"

        [ ! -s "${FILTERED}" ] && continue

        PEAK_COUNT=$(wc -l < "${FILTERED}")

        computeMatrix reference-point \
            -S ${BW_FILES} -R "${FILTERED}" \
            --referencePoint center -a 3000 -b 3000 \
            --skipZeros --missingDataAsZero \
            -o "${GC_DIR}/${marker}_peaks_matrix.gz" \
            --samplesLabel ${LABEL_ARRAY[@]} -p ${THREADS}

        plotHeatmap -m "${GC_DIR}/${marker}_peaks_matrix.gz" \
            -out "${GC_DIR}/${marker}_${gc}_peaks_heatmap.pdf" \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 15 --heatmapWidth 10 \
            --plotTitle "${marker} at ${gc}-Associated Peaks (n=${PEAK_COUNT})" \
            --refPointLabel "Peak center" --sortUsing mean --sortRegions descend
    done
done

echo "Gene-specific analysis complete! Results in: ${OUT_DIR}"
