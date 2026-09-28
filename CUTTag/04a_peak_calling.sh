#!/bin/bash
#SBATCH --job-name=CUTTag_PeakCalling
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 04a: Peak Calling with MACS2 (Treatment vs Untreated control)
# Calls peaks for each condition using untreated samples as control
#
# Input:  Deduplicated BAM files from Step 03
# Output: Peak files (narrowPeak format), peak summaries
#==================================================================================

source $(dirname "$0")/00_config.sh

echo "Starting MACS2 peak calling..."

# Process each replicate separately
for rep in ${REPLICATES[@]}; do
    for condition in ${CONDITIONS[@]}; do
        for marker in ${MARKERS[@]}; do
            echo "Processing: ${condition} ${marker} rep${rep} vs ${TIME_CONTROL}"

            # Create output directory
            outdir="${PEAK_DIR}/${ORGANISM}/${TIME_CONTROL}_control/${condition}/${marker}/rep${rep}"
            mkdir -p ${outdir}

            # Get BAM files from sample map
            treatment_bam="${DEDUP_DIR}/${sample_map[${condition}_${marker}_${rep}]}_dedup.bam"
            control_bam="${DEDUP_DIR}/${sample_map[${TIME_CONTROL}_${marker}_${rep}]}_dedup.bam"

            if [ -f "$treatment_bam" ] && [ -f "$control_bam" ]; then
                macs2 callpeak \
                    -t ${treatment_bam} \
                    -c ${control_bam} \
                    -g ${GENOME_CODE} \
                    -f ${PEAK_FORMAT} \
                    -n macs2_peak_q${MACS2_QVALUE} \
                    --outdir ${outdir} \
                    -q ${MACS2_QVALUE} \
                    --keep-dup all \
                    2>${outdir}/macs2Peak_summary.txt

                # Rename narrowPeak to .bed
                peak_file="${outdir}/macs2_peak_q${MACS2_QVALUE}_peaks.narrowPeak"
                if [ -f "$peak_file" ]; then
                    mv ${peak_file} ${peak_file}.bed
                fi
            else
                echo "  WARNING: Missing BAM files for ${condition} ${marker} rep${rep}"
            fi
        done
    done
done

# Generate summary report
summary_file="${PEAK_DIR}/peak_summary.txt"
echo -e "Replicate\tComparison\tTarget\tQ_value\tMode\tPeak_Count" > ${summary_file}

for rep in ${REPLICATES[@]}; do
    for condition in ${CONDITIONS[@]}; do
        for marker in ${MARKERS[@]}; do
            peak_file="${PEAK_DIR}/${ORGANISM}/${TIME_CONTROL}_control/${condition}/${marker}/rep${rep}/macs2_peak_q${MACS2_QVALUE}_peaks.narrowPeak.bed"
            if [ -f "$peak_file" ]; then
                peak_count=$(wc -l < "$peak_file")
                echo -e "Rep${rep}\t${condition}_vs_${TIME_CONTROL}\t${marker}\t${MACS2_QVALUE}\tdefault\t${peak_count}" >> ${summary_file}
            fi
        done
    done
done

echo "Peak calling complete! Summary: ${summary_file}"
