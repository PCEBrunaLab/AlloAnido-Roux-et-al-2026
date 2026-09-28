#!/bin/bash
#SBATCH --job-name=CUTTag_QC
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 02: Quality Control
# Generates QC metrics (alignment rate, duplication, library complexity)
#
# Input:  Sorted BAM files from Step 01
# Output: QC summary table, QC plots (PDF)
#==================================================================================

source $(dirname "$0")/00_config.sh

mkdir -p ${QC_DIR} ${LOG_DIR}

# Create summary header
echo -e "Sample\tTotal_Reads\tMapped_Reads\tAlignment_Rate\tDuplication_Rate\tUnique_Fragments\tSeq_Depth_M\tMapped_Frags_M\tCondition\tTarget" > ${QC_DIR}/qc_summary.txt

# Process each sample
while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
    SAMPLE=$(echo "$SAMPLE" | tr -d '[:space:]')

    ORIGINAL_BAM=${BAM_DIR}/${SAMPLE}_sorted.bam

    if [ ! -f "$ORIGINAL_BAM" ]; then
        echo "BAM not found for ${SAMPLE}, skipping..."
        continue
    fi

    # Calculate metrics
    TOTAL_READS=$(samtools view -c ${ORIGINAL_BAM})
    MAPPED_READS=$(samtools view -c -F 4 ${ORIGINAL_BAM})
    ALIGNMENT_RATE=$(echo "scale=1; (${MAPPED_READS} * 100)/${TOTAL_READS}" | bc)
    SEQ_DEPTH=$(echo "scale=2; ${TOTAL_READS}/1000000" | bc)
    MAPPED_FRAGS=$(echo "scale=2; (${MAPPED_READS}/2)/1000000" | bc)

    # Duplication metrics (Picard)
    picard MarkDuplicates \
        I=${ORIGINAL_BAM} \
        O=${QC_DIR}/${SAMPLE}_marked.bam \
        M=${QC_DIR}/${SAMPLE}_dup_metrics.txt \
        REMOVE_DUPLICATES=false \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR=/tmp \
        QUIET=true

    if [ -f "${QC_DIR}/${SAMPLE}_dup_metrics.txt" ]; then
        DUP_RATE=$(grep -A 1 "LIBRARY" ${QC_DIR}/${SAMPLE}_dup_metrics.txt | tail -1 | cut -f9 | awk '{printf "%.1f", $1*100}')
        [ -z "$DUP_RATE" ] && DUP_RATE="0"
    else
        DUP_RATE="0"
    fi

    # Unique fragments
    if [ -f "${QC_DIR}/${SAMPLE}_marked.bam" ]; then
        UNIQUE_FRAGS=$(samtools view -c -f 2 -F 1024 ${QC_DIR}/${SAMPLE}_marked.bam | awk '{print int($1/2)}')
    else
        UNIQUE_FRAGS="0"
    fi

    # Determine condition and target from sample metadata
    # NOTE: Adapt this section to your sample naming convention
    CONDITION="Unknown"
    TARGET="Unknown"
    # Example: parse from sample_list metadata or use a lookup table

    echo -e "${SAMPLE}\t${TOTAL_READS}\t${MAPPED_READS}\t${ALIGNMENT_RATE}\t${DUP_RATE}\t${UNIQUE_FRAGS}\t${SEQ_DEPTH}\t${MAPPED_FRAGS}\t${CONDITION}\t${TARGET}" >> ${QC_DIR}/qc_summary.txt

    # Clean up temporary files
    rm -f ${QC_DIR}/${SAMPLE}_marked.bam ${QC_DIR}/${SAMPLE}_marked.bam.bai

done < "${SAMPLE_LIST}"

# Generate QC plots with R
cat > ${QC_DIR}/plot_qc.R << 'RSCRIPT'
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

data <- read.table("qc_summary.txt", header=TRUE, sep="\t")

target_colors <- c("H3K4me3"="#E41A1C", "IgG"="#377EB8", "KDM5A"="#4DAF4A", "KDM5B"="#984EA3")
data$Condition.Mark <- interaction(data$Condition, data$Target, sep=".")

p1 <- ggplot(data, aes(x=Condition.Mark, y=Seq_Depth_M, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    labs(title="Sequencing Depth", y="Sequencing Depth (M)", x="") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))

p2 <- ggplot(data, aes(x=Condition.Mark, y=Duplication_Rate, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    geom_hline(yintercept=40, linetype="dashed", color="orange", alpha=0.5) +
    geom_hline(yintercept=60, linetype="dashed", color="red", alpha=0.5) +
    labs(title="Duplication Rate", y="Duplication Rate (%)", x="") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))

p3 <- ggplot(data, aes(x=Condition.Mark, y=Unique_Fragments/1000, fill=Target)) +
    geom_boxplot(alpha=0.8, width=0.7) +
    geom_point(position=position_jitter(width=0.15), size=1.5, alpha=0.8) +
    scale_fill_manual(values=target_colors, name="Mark") +
    labs(title="Unique Library Size", y="Unique Fragments (K)", x="") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))

pdf("CUT_Tag_QC_Combined.pdf", width=12, height=12)
grid.arrange(p1, p2, p3, ncol=1)
dev.off()

cat("QC plots saved successfully!\n")
RSCRIPT

cd ${QC_DIR} && Rscript plot_qc.R

echo "QC analysis complete! Results in: ${QC_DIR}"
