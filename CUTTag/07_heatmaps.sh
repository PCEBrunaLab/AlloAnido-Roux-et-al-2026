#!/bin/bash
#SBATCH --job-name=CUTTag_Heatmaps
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 07: Genome-wide Heatmaps (TSS and Peak-centered)
# Creates time-series heatmaps with deepTools (computeMatrix + plotHeatmap)
#
# Input:  BigWig files from Step 06, gene annotations, merged peak files
# Output: Heatmap PDFs, profile plots, deepTools matrices
#==================================================================================

source $(dirname "$0")/00_config.sh

OUT_DIR="${HEATMAP_DIR}/final_timeseries"
mkdir -p ${OUT_DIR}

# Create TSS file from genes.bed
TSS_FILE="${PROJ_DIR}/genes_tss.bed"
if [ -f "${GENES_BED}" ]; then
    awk 'BEGIN{OFS="\t"} {
        if(NF < 6) next;
        if($6=="+") start=$2;
        else if($6=="-") start=$3-1;
        else next;
        if(start<0) start=0;
        print $1,start,start+1,$4"_TSS",0,$6
    }' ${GENES_BED} | sort -k1,1 -k2,2n | uniq > ${TSS_FILE}
    echo "TSS file created: $(wc -l < ${TSS_FILE}) regions"
else
    echo "ERROR: genes.bed not found at ${GENES_BED}"
    exit 1
fi

# Process each marker and merge distance
for marker in ${MARKERS[@]}; do
    echo "Processing ${marker}..."

    for distance in ${HOMER_MERGE_DISTANCES[@]}; do
        PEAK_FILE="${PROJ_DIR}/Homer/mergePeaks/${ORGANISM}/${TIME_CONTROL}_control/${marker}/${distance}/mergedPeakFile.bed"

        [ ! -f "$PEAK_FILE" ] && continue
        PEAK_COUNT=$(wc -l < "$PEAK_FILE")

        # Build BigWig file list
        BW_FILES=""
        for condition in ${ALL_CONDITIONS[@]}; do
            bw="${BIGWIG_DIR}/${ORGANISM}/${condition}/${marker}/normalized.bw"
            [ -f "$bw" ] && BW_FILES="${BW_FILES} ${bw}"
        done

        [ -z "$BW_FILES" ] && continue

        # TSS-centered heatmap
        computeMatrix reference-point \
            -S ${BW_FILES} \
            -R ${TSS_FILE} \
            --referencePoint TSS \
            -a ${HEATMAP_DISTANCE} -b ${HEATMAP_DISTANCE} \
            --skipZeros --missingDataAsZero \
            -o ${OUT_DIR}/${marker}_TSS_merged_${distance}bp_matrix.gz \
            --samplesLabel ${CONDITION_LABELS[@]} \
            -p ${THREADS}

        plotHeatmap \
            -m ${OUT_DIR}/${marker}_TSS_merged_${distance}bp_matrix.gz \
            -out ${OUT_DIR}/${marker}_TSS_merged_${distance}bp.pdf \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 12 --heatmapWidth 8 \
            --plotTitle "${marker} - TSS Signal (${distance}bp Merged)" \
            --xAxisLabel "Distance from TSS (bp)" \
            --refPointLabel "TSS" \
            --sortUsing mean --sortRegions descend \
            --zMin 0 --zMax auto

        # Peak-centered heatmap
        computeMatrix reference-point \
            -S ${BW_FILES} \
            -R ${PEAK_FILE} \
            --referencePoint center \
            -a ${HEATMAP_DISTANCE} -b ${HEATMAP_DISTANCE} \
            --skipZeros --missingDataAsZero \
            -o ${OUT_DIR}/${marker}_Peak_merged_${distance}bp_matrix.gz \
            --samplesLabel ${CONDITION_LABELS[@]} \
            -p ${THREADS}

        plotHeatmap \
            -m ${OUT_DIR}/${marker}_Peak_merged_${distance}bp_matrix.gz \
            -out ${OUT_DIR}/${marker}_Peak_merged_${distance}bp.pdf \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 12 --heatmapWidth 8 \
            --plotTitle "${marker} - Peak Signal (${distance}bp Merged)" \
            --xAxisLabel "Distance from Peak Center (bp)" \
            --yAxisLabel "Peaks (n=${PEAK_COUNT})" \
            --refPointLabel "Center" \
            --sortUsing mean --sortRegions descend \
            --zMin 0 --zMax auto
    done
done

echo "Heatmap generation complete! Results in: ${OUT_DIR}"
