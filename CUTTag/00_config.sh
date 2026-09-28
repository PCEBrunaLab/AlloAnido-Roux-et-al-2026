#!/bin/bash
#==================================================================================
# CUT&Tag Pipeline Configuration
# Edit this file to set paths and parameters for your environment
#==================================================================================

# Project directory (set to your working directory)
export PROJ_DIR="/path/to/your/project"

# Input data
export FASTQ_DIR="${PROJ_DIR}/raw_fastq"
export SAMPLE_LIST="${PROJ_DIR}/sample_list.txt"

# Reference genome
export GENOME_DIR="/path/to/genome"
export BOWTIE2_INDEX="${GENOME_DIR}/bowtie2_index/GRCh38"
export GENOME_SIZE="2913022398"  # Human GRCh38 effective genome size for RPGC
export GENOME_CODE="hs"           # MACS2 genome code (hs=human, mm=mouse)

# Output directories (created automatically)
export BAM_DIR="${PROJ_DIR}/bam_files"
export DEDUP_DIR="${PROJ_DIR}/bam_files_dedup"
export QC_DIR="${PROJ_DIR}/QC_final"
export PEAK_DIR="${PROJ_DIR}/peakCalling/MACS2"
export BIGWIG_DIR="${PROJ_DIR}/alignment/bigwig"
export HEATMAP_DIR="${PROJ_DIR}/heatmap"
export ANALYSIS_DIR="${PROJ_DIR}/Analysis"
export LOG_DIR="${PROJ_DIR}/logs"

# Reference files
export GENES_BED="${PROJ_DIR}/genes.bed"
export GENE_LIST_CSV="${PROJ_DIR}/Reference/gene_list.csv"

# Experimental design
export ORGANISM="human"
export TIME_CONTROL="Untreated"
export CONDITIONS=("D7_Cis" "D14_Cis" "D14_C70")
export MARKERS=("H3K4me3" "KDM5A" "KDM5B")
export REPLICATES=(1 2)

# Condition labels for plots
export CONDITION_LABELS=("Day 0" "Day 7" "Day 14" "Day 14+C70")
export ALL_CONDITIONS=("Untreated" "D7_Cis" "D14_Cis" "D14_C70")

# SLURM settings (adjust for your cluster)
export SLURM_PARTITION="compute"
export SLURM_MAIL="your.email@institution.edu"

# Peak calling parameters
export MACS2_QVALUE=0.1
export PEAK_FORMAT="BAMPE"

# BigWig parameters
export BIGWIG_BINSIZE=10
export BIGWIG_NORM="RPGC"

# Heatmap parameters
export HEATMAP_DISTANCE=3000  # bp upstream/downstream of reference point
export HOMER_MERGE_DISTANCES=(5000 7500 10000)

# Peak-to-gene assignment
export PEAK_GENE_DISTANCE=5000  # Max distance for peak-to-gene assignment

# Number of threads
export THREADS=8

#==================================================================================
# Sample mapping - EDIT THIS FOR YOUR SAMPLES
# Format: sample_map["Condition_Marker_Replicate"]="SampleID"
#==================================================================================
declare -A sample_map

# Example mapping (replace with your sample IDs):
# Replicate 1
# sample_map["Untreated_H3K4me3_1"]="Sample_001"
# sample_map["D7_Cis_H3K4me3_1"]="Sample_002"
# ...

# Replicate 2
# sample_map["Untreated_H3K4me3_2"]="Sample_007"
# sample_map["D7_Cis_H3K4me3_2"]="Sample_008"
# ...

export sample_map
