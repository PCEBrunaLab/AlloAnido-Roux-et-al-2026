#!/bin/bash
#SBATCH --job-name=CUTTag_CellState
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 09: Cell State Gene Signature Analysis
# Analyzes CUT&Tag signal at ADRN/MES/Intermediate cell state gene signatures
# Creates TSS, gene body, and peak-centered heatmaps per cell state
#
# Input:  BigWig files (Step 06), cell state gene lists, gene annotations
# Output: Cell state-specific heatmaps, active gene analysis
#==================================================================================

source $(dirname "$0")/00_config.sh

OUT_DIR="${HEATMAP_DIR}/cellstate_genelist"
mkdir -p ${OUT_DIR}

HOMER_DIST=7500

# Cell state gene files (edit paths for your data)
ADRN_FILE="${PROJ_DIR}/Reference/vanGroningen_2017_ADRN.csv"
MES_FILE="${PROJ_DIR}/Reference/vanGroningen_2017_MES.csv"
INT_FILE="${PROJ_DIR}/Reference/Int.csv"

echo "Cell State Gene-Specific Analysis"

# Parse cell state gene lists
declare -A GENE_LISTS

for state_file in "${ADRN_FILE}:ADRN" "${MES_FILE}:MES" "${INT_FILE}:Int"; do
    IFS=':' read -r filepath state_name <<< "$state_file"
    [ ! -f "$filepath" ] && continue

    OUT_LIST="${OUT_DIR}/${state_name}_genes.txt"
    if head -1 "$filepath" | grep -q "Gene"; then
        tail -n +2 "$filepath" | awk -F',' '{print $1}' | tr -d ' \r' | sort -u > "$OUT_LIST"
    else
        awk -F',' '{print $1}' "$filepath" | tr -d ' \r' | sort -u > "$OUT_LIST"
    fi

    count=$(wc -l < "$OUT_LIST")
    echo "  ${state_name}: ${count} genes"
    GENE_LISTS["${state_name}"]="$OUT_LIST"
done

# Process each cell state
for state in "${!GENE_LISTS[@]}"; do
    STATE_DIR="${OUT_DIR}/${state}"
    mkdir -p "${STATE_DIR}"

    GENE_LIST="${GENE_LISTS[${state}]}"

    SUBSET_BED="${STATE_DIR}/genes_subset.bed"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${GENE_LIST}" "${GENES_BED}" | sort -k1,1 -k2,2n | uniq > "${SUBSET_BED}"

    NUM_GENES=$(wc -l < "${SUBSET_BED}" 2>/dev/null || echo "0")
    [ "$NUM_GENES" -eq 0 ] && continue

    # TSS BED
    TSS_BED="${STATE_DIR}/tss.bed"
    awk 'BEGIN{OFS="\t"} {
        if(NF<6) next;
        if($6=="+") tss=$2; else if($6=="-") tss=$3-1; else next;
        if(tss<0) tss=0;
        print $1, tss, tss+1, $4"_TSS", 0, $6
    }' "${SUBSET_BED}" | sort -k1,1 -k2,2n | uniq > "${TSS_BED}"

    for marker in ${MARKERS[@]}; do
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
            -o "${STATE_DIR}/${marker}_TSS_matrix.gz" \
            --samplesLabel ${LABEL_ARRAY[@]} -p ${THREADS}

        plotHeatmap -m "${STATE_DIR}/${marker}_TSS_matrix.gz" \
            -out "${STATE_DIR}/${marker}_${state}_TSS_heatmap.pdf" \
            --colorMap RdBu_r --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 15 --heatmapWidth 10 \
            --plotTitle "${marker} at ${state} Cell State Gene TSS" \
            --refPointLabel "TSS" --sortUsing mean --sortRegions descend

        # Gene body (scale-regions)
        computeMatrix scale-regions \
            -S ${BW_FILES} -R "${SUBSET_BED}" \
            --beforeRegionStartLength 3000 \
            --regionBodyLength 5000 \
            --afterRegionStartLength 3000 \
            --skipZeros --missingDataAsZero \
            -o "${STATE_DIR}/${marker}_gene_body_matrix.gz" \
            --samplesLabel ${LABEL_ARRAY[@]} -p ${THREADS}

        plotHeatmap -m "${STATE_DIR}/${marker}_gene_body_matrix.gz" \
            -out "${STATE_DIR}/${marker}_${state}_gene_body_heatmap.pdf" \
            --colorMap RdBu_r --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 15 --heatmapWidth 10 \
            --plotTitle "${marker} Across ${state} Gene Bodies" \
            --startLabel "Start" --endLabel "End" \
            --sortUsing mean --sortRegions descend
    done
done

echo "Cell state analysis complete! Results in: ${OUT_DIR}"
