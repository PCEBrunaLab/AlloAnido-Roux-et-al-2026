#!/bin/bash
#SBATCH --job-name=CUTTag_BigWig
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 06: BigWig Generation (RPGC normalized)
# Creates normalized coverage tracks for genome browser visualization
#
# Input:  Deduplicated BAM files from Step 03
# Output: RPGC-normalized BigWig files (.bw) per replicate and merged
#==================================================================================

source $(dirname "$0")/00_config.sh

echo "Creating RPGC Normalized BigWig Files..."

for condition in ${ALL_CONDITIONS[@]}; do
    for marker in ${MARKERS[@]}; do
        # Per-replicate BigWigs
        for rep in ${REPLICATES[@]}; do
            sample_key="${condition}_${marker}_${rep}"
            sample_id="${sample_map[$sample_key]}"

            [ -z "$sample_id" ] && continue

            mkdir -p ${BIGWIG_DIR}/${ORGANISM}/${condition}/${marker}

            input_bam="${DEDUP_DIR}/${sample_id}_dedup.bam"
            output_bw="${BIGWIG_DIR}/${ORGANISM}/${condition}/${marker}/normalized_rep${rep}.bw"

            [ ! -f "$input_bam" ] && continue

            echo "Processing: ${condition} ${marker} rep${rep}"

            bamCoverage --normalizeUsing ${BIGWIG_NORM} \
                --effectiveGenomeSize ${GENOME_SIZE} \
                -b ${input_bam} \
                -o ${output_bw} \
                --binSize ${BIGWIG_BINSIZE} \
                --numberOfProcessors ${THREADS}
        done

        # Merged BigWig (both replicates)
        bam1="${DEDUP_DIR}/${sample_map[${condition}_${marker}_1]}_dedup.bam"
        bam2="${DEDUP_DIR}/${sample_map[${condition}_${marker}_2]}_dedup.bam"

        if [ -f "$bam1" ] && [ -f "$bam2" ]; then
            merged_bam="${BIGWIG_DIR}/${ORGANISM}/${condition}/${marker}/merged.bam"
            samtools merge -@ ${THREADS} ${merged_bam} ${bam1} ${bam2}
            samtools index ${merged_bam}

            bamCoverage --normalizeUsing ${BIGWIG_NORM} \
                --effectiveGenomeSize ${GENOME_SIZE} \
                -b ${merged_bam} \
                -o ${BIGWIG_DIR}/${ORGANISM}/${condition}/${marker}/normalized.bw \
                --binSize ${BIGWIG_BINSIZE} \
                --numberOfProcessors ${THREADS}

            rm -f ${merged_bam} ${merged_bam}.bai
        fi
    done
done

echo "BigWig generation complete!"
