#!/bin/bash
#SBATCH --job-name=CUTTag_Mapping
#SBATCH --partition=${SLURM_PARTITION:-compute}
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --array=1-48%8
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

#==================================================================================
# Step 01: Read Mapping (Bowtie2)
# Maps paired-end CUT&Tag FASTQ files to reference genome
#
# Input:  Paired-end FASTQ files (.fastq.gz)
# Output: Sorted BAM files (*_sorted.bam), alignment statistics
#==================================================================================

# Source configuration
source $(dirname "$0")/00_config.sh

# Create output directories
mkdir -p ${BAM_DIR} ${LOG_DIR} ${BAM_DIR}/qc

# Get the specific sample for this array job
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample found for array index ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Define input files
R1=$(ls ${FASTQ_DIR}/${SAMPLE}_S*_L*_R1_001.fastq.gz 2>/dev/null | head -n1)
R2=$(ls ${FASTQ_DIR}/${SAMPLE}_S*_L*_R2_001.fastq.gz 2>/dev/null | head -n1)

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Error: Input files not found for sample ${SAMPLE}"
    exit 1
fi

# Skip if output already exists
if [ -f "${BAM_DIR}/${SAMPLE}_sorted.bam" ]; then
    echo "[$(date)] Skipping ${SAMPLE} - output already exists"
    exit 0
fi

echo "Processing sample: ${SAMPLE}"

# Map with CUT&Tag-specific parameters
bowtie2 --local \
    --very-sensitive \
    --no-unal \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    -p ${THREADS} \
    -x ${BOWTIE2_INDEX} \
    -1 ${R1} \
    -2 ${R2} \
    2> ${LOG_DIR}/${SAMPLE}.bowtie2.log | \
    samtools sort -@ ${THREADS} -O bam -o ${BAM_DIR}/${SAMPLE}_sorted.bam

if [ $? -ne 0 ]; then
    echo "Error: Mapping failed for ${SAMPLE}"
    exit 1
fi

# Index BAM
samtools index -@ ${THREADS} ${BAM_DIR}/${SAMPLE}_sorted.bam

# Alignment statistics
samtools flagstat ${BAM_DIR}/${SAMPLE}_sorted.bam > ${BAM_DIR}/qc/${SAMPLE}.flagstat.txt

echo "[$(date)] Completed mapping for: ${SAMPLE}"
