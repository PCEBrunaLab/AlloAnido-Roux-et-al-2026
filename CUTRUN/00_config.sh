#!/bin/bash
# ============================================================
# Title:        Shared configuration for the SK-N-SH CUT&RUN analysis
# Input:        samplesheet.csv, sample_list.txt
# Output:       exported variables, sample maps and helper functions
# Depends on:   bash >= 4
# Notes:        Sourced by every script in this directory. Edit PROJ_DIR,
#               GENOME_DIR and CUTTAG_DIR before the first run.
#               Sample names follow Condition_Target_Replicate
#               (for example Awakening_H3K4me3_1).
# ============================================================

# ---- 1. Paths --------------------------------------------------------------

export PROJ_DIR="${PROJ_DIR:-/path/to/project}"
export GENOME_DIR="${GENOME_DIR:-/path/to/genome}"
export CUTTAG_DIR="${CUTTAG_DIR:-${PROJ_DIR}/CUTTag}"      # parental SK-N-SH CUT&Tag (BAM_dedup/)

export FASTQ_RUN_DIR="${PROJ_DIR}/raw_fastq"               # <sample>_R{1,2}_001.fastq.gz
export TRIM_DIR="${PROJ_DIR}/fastq_trimmed"
export SAMPLESHEET="${PROJ_DIR}/samplesheet.csv"
export SAMPLE_LIST="${PROJ_DIR}/sample_list.txt"

export BOWTIE2_INDEX="${GENOME_DIR}/bowtie2_index/GRCh38"
export REF_DIR="${PROJ_DIR}/reference"
export CHROM_SIZES="${REF_DIR}/hg38.chrom.sizes"

export BAM_DIR="${PROJ_DIR}/bam_files"
export DEDUP_DIR="${PROJ_DIR}/bam_files_dedup"
export MERGED_BAM_DIR="${PROJ_DIR}/merged_bam_files"
export DEDUP_QC_DIR="${PROJ_DIR}/qc"
export RESULTS_DIR="${PROJ_DIR}/results"
export LOG_DIR="${PROJ_DIR}/logs"

# ---- 2. Design -------------------------------------------------------------

export ALL_CONDITIONS=("Awakening" "C70" "T5")
export CONDITION_LABELS=("Awakening" "C70" "T5")
export MARKERS=("H3K4me3" "KDM5A" "ATF2" "Fos" "Jun")
export REPLICATES=(1 2)

# ---- 3. Parameters ---------------------------------------------------------

export THREADS=8
export MIN_READ_LENGTH=20          # Cutadapt -m
export GENOME_CODE="hs"            # MACS3 -g
export PEAK_FORMAT="BAMPE"         # MACS3 -f
export MACS_QVALUE=0.1             # MACS3 -q
export HEATMAP_BINSIZE=50
export HEATMAP_DISTANCE=3000

# Optional named environments.
export MACS_ENV="${MACS_ENV:-}"
export DEEPTOOLS_ENV="${DEEPTOOLS_ENV:-}"
export DEEPTOOLS_BIN="${DEEPTOOLS_BIN:-}"

# ENCODE hg38 blacklist v2 (Amemiya, Kundaje & Boyle 2019).
export BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
export BLACKLIST_BED="${REF_DIR}/hg38-blacklist.v2.bed"
export BLACKLIST_BED3="${REF_DIR}/hg38-blacklist.v2.bed3"

# ---- 4. Sample map ---------------------------------------------------------
# samplesheet.csv columns: sample_id,sample_name,condition,target,replicate

declare -A sample_map sample_condition sample_target sample_replicate

if [ -f "${SAMPLESHEET}" ]; then
    while IFS=, read -r sid sname cond targ rep _rest; do
        [ "${sid}" = "sample_id" ] && continue
        [ -z "${sname}" ] && continue
        sample_map["${cond}_${targ}_${rep}"]="${sname}"
        sample_condition["${sname}"]="${cond}"
        sample_target["${sname}"]="${targ}"
        sample_replicate["${sname}"]="${rep}"
    done < "${SAMPLESHEET}"
else
    echo "WARNING: samplesheet not found at ${SAMPLESHEET}" >&2
fi

# ---- 5. Helpers ------------------------------------------------------------

# Put deepTools on PATH (module, DEEPTOOLS_BIN or DEEPTOOLS_ENV) and check it runs.
activate_deeptools () {
    local tool="${1:-computeMatrix}"
    if ! command -v "${tool}" >/dev/null 2>&1; then
        if type module >/dev/null 2>&1; then
            module load deeptools 2>/dev/null || module load deepTools 2>/dev/null || true
        fi
        [ -n "${DEEPTOOLS_BIN}" ] && export PATH="${DEEPTOOLS_BIN}:${PATH}"
        if [ -n "${DEEPTOOLS_ENV}" ] && command -v conda >/dev/null 2>&1; then
            set +u; conda activate "${DEEPTOOLS_ENV}" 2>/dev/null || true; set -u
        fi
    fi
    "${tool}" --version >/dev/null 2>&1 || { echo "ERROR: deepTools '${tool}' not available" >&2; return 1; }
    echo "deepTools: $(${tool} --version 2>&1 | tr -d '\n')"
}

# Put MACS3 on PATH (module or MACS_ENV) and print its path.
find_macs () {
    if type module >/dev/null 2>&1; then
        module load macs3 2>/dev/null || module load MACS3 2>/dev/null || true
    fi
    if [ -n "${MACS_ENV}" ] && command -v conda >/dev/null 2>&1; then
        set +u; conda activate "${MACS_ENV}" 2>/dev/null || true; set -u
    fi
    command -v macs3 >/dev/null 2>&1 && macs3 --version >/dev/null 2>&1 \
        || { echo "ERROR: macs3 not available" >&2; return 1; }
    command -v macs3
}

# Return the chromosome sizes file, writing it from a BAM header if missing.
ensure_chrom_sizes () {
    if [ -s "${CHROM_SIZES}" ]; then echo "${CHROM_SIZES}"; return 0; fi
    local bam
    bam=$(ls "${DEDUP_DIR}"/*_dedup.bam 2>/dev/null | head -1)
    [ -n "${bam}" ] || { echo "ERROR: no chrom.sizes and no BAM to derive them" >&2; return 1; }
    mkdir -p "${REF_DIR}"
    samtools view -H "${bam}" | awk '/^@SQ/{sn="";ln="";
        for(i=2;i<=NF;i++){ if($i~/^SN:/) sn=substr($i,4); if($i~/^LN:/) ln=substr($i,4) }
        if(sn!="" && ln!="") print sn"\t"ln }' > "${CHROM_SIZES}"
    [ -s "${CHROM_SIZES}" ] && echo "${CHROM_SIZES}"
}

# Rename the chromosomes of a BED file to the style used in a BAM ("1" or "chr1").
match_chrom_to_bam () {
    local src="$1" bam="$2" dest="$3" bam_style src_style
    bam_style=$(samtools view -H "${bam}" | awk '/^@SQ/{sub(/^SN:/,"",$2); print $2; exit}')
    src_style=$(awk 'NF && $1 !~ /^#/ {print $1; exit}' "${src}")
    if [[ "${bam_style}" == chr* && "${src_style}" != chr* ]]; then
        awk 'BEGIN{OFS="\t"} NF && $1 !~ /^#/ {$1="chr"$1; print}' "${src}" > "${dest}"
    elif [[ "${bam_style}" != chr* && "${src_style}" == chr* ]]; then
        awk 'BEGIN{OFS="\t"} NF && $1 !~ /^#/ {sub(/^chr/,"",$1); print}' "${src}" > "${dest}"
    else
        awk 'NF && $1 !~ /^#/' "${src}" > "${dest}"
    fi
}

mkdir -p "${TRIM_DIR}" "${BAM_DIR}" "${DEDUP_DIR}" "${MERGED_BAM_DIR}" "${DEDUP_QC_DIR}" \
         "${RESULTS_DIR}" "${REF_DIR}" "${LOG_DIR}" 2>/dev/null || true
