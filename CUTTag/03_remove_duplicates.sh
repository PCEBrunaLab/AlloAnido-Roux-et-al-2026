#!/bin/bash
#SBATCH --job-name=CUTTag_RemoveDup
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 03: Duplicate Removal (Picard MarkDuplicates)
# Removes PCR duplicates and generates before/after comparison metrics
#
# Input:  Sorted BAM files from Step 01
# Output: Deduplicated BAM files (*_dedup.bam), deduplication metrics
#==================================================================================

source $(dirname "$0")/00_config.sh

METRICS_DIR=${DEDUP_DIR}/metrics
COMP_DIR=${PROJ_DIR}/QC_dedup_comparison

mkdir -p ${DEDUP_DIR} ${METRICS_DIR} ${LOG_DIR} ${COMP_DIR}

# Initialize post-dedup summary
echo -e "Sample\tTotal_Reads_PostDedup\tMapped_PostDedup\tDuplicates_Removed\tPercent_Dup_Removed\tUnique_Fragments_Final\tSeq_Depth_PostDedup_M\tCondition\tTarget" > ${COMP_DIR}/qc_summary_postdedup.txt

while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
    SAMPLE=$(echo "$SAMPLE" | tr -d '[:space:]')

    INPUT_BAM=${BAM_DIR}/${SAMPLE}_sorted.bam
    OUTPUT_BAM=${DEDUP_DIR}/${SAMPLE}_dedup.bam
    METRICS_FILE=${METRICS_DIR}/${SAMPLE}.dedup_metrics.txt

    if [ ! -f "$INPUT_BAM" ]; then
        echo "Input BAM not found for ${SAMPLE}, skipping..."
        continue
    fi

    if [ -f "$OUTPUT_BAM" ]; then
        echo "Output exists for ${SAMPLE}, skipping deduplication..."
    else
        echo "Removing duplicates from ${SAMPLE}..."

        picard MarkDuplicates \
            I=${INPUT_BAM} \
            O=${OUTPUT_BAM} \
            M=${METRICS_FILE} \
            REMOVE_DUPLICATES=true \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=/tmp

        if [ $? -ne 0 ]; then
            echo "Error: Deduplication failed for ${SAMPLE}"
            continue
        fi

        samtools index ${OUTPUT_BAM}
    fi

    # Post-dedup metrics
    TOTAL_POST=$(samtools view -c ${OUTPUT_BAM})
    MAPPED_POST=$(samtools view -c -F 4 ${OUTPUT_BAM})
    SEQ_DEPTH_POST=$(echo "scale=2; ${TOTAL_POST}/1000000" | bc)
    UNIQUE_FINAL=$(samtools view -c -f 2 ${OUTPUT_BAM} | awk '{print int($1/2)}')

    if [ -f "$METRICS_FILE" ]; then
        DUPS_REMOVED=$(grep -A 1 "LIBRARY" ${METRICS_FILE} | tail -1 | cut -f 8)
        PERCENT_DUP=$(grep -A 1 "LIBRARY" ${METRICS_FILE} | tail -1 | cut -f 9 | awk '{printf "%.1f", $1*100}')
    else
        DUPS_REMOVED="NA"
        PERCENT_DUP="NA"
    fi

    # NOTE: Adapt condition/target parsing to your sample naming
    CONDITION="Unknown"
    TARGET="Unknown"

    echo -e "${SAMPLE}\t${TOTAL_POST}\t${MAPPED_POST}\t${DUPS_REMOVED}\t${PERCENT_DUP}\t${UNIQUE_FINAL}\t${SEQ_DEPTH_POST}\t${CONDITION}\t${TARGET}" >> ${COMP_DIR}/qc_summary_postdedup.txt

done < "${SAMPLE_LIST}"

echo "Deduplication complete! Results in: ${DEDUP_DIR}"
