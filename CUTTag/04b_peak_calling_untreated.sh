#!/bin/bash
#SBATCH --job-name=CUTTag_PeakCalling_Untreated
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

#==================================================================================
# Step 04b: Peak Calling for Untreated samples (vs IgG control)
# Creates baseline peaks for the untreated condition
#
# Input:  Deduplicated BAM files from Step 03
# Output: Untreated peak files (narrowPeak format)
#==================================================================================

source $(dirname "$0")/00_config.sh

outBase="${PEAK_DIR}/${ORGANISM}/IgG_control/${TIME_CONTROL}"

for marker in ${MARKERS[@]}; do
    echo "--- Processing ${marker} ---"

    for rep in ${REPLICATES[@]}; do
        out_dir="${outBase}/${marker}/rep${rep}"
        mkdir -p ${out_dir}

        treatment_bam="${DEDUP_DIR}/${sample_map[${TIME_CONTROL}_${marker}_${rep}]}_dedup.bam"
        control_bam="${DEDUP_DIR}/${sample_map[${TIME_CONTROL}_IgG_${rep}]}_dedup.bam"

        if [ ! -f "$treatment_bam" ] || [ ! -f "$control_bam" ]; then
            echo "  Missing BAM files for ${marker} rep${rep}, skipping..."
            continue
        fi

        macs2 callpeak \
            -t ${treatment_bam} \
            -c ${control_bam} \
            -g ${GENOME_CODE} \
            -f ${PEAK_FORMAT} \
            -n macs2_peak_q${MACS2_QVALUE} \
            --outdir ${out_dir} \
            -q ${MACS2_QVALUE} \
            --keep-dup all \
            2>${out_dir}/macs2Peak_summary.txt

        # Rename to .bed
        peak_file="${out_dir}/macs2_peak_q${MACS2_QVALUE}_peaks.narrowPeak"
        if [ -f "$peak_file" ]; then
            mv ${peak_file} ${peak_file}.bed
            peak_count=$(wc -l < ${peak_file}.bed)
            echo "  ${marker} rep${rep}: ${peak_count} peaks"
        fi
    done
done

echo "Untreated peak calling complete! Output: ${outBase}"
