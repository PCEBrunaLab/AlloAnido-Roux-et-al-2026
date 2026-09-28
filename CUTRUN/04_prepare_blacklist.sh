#!/bin/bash
# ============================================================
# Title:        Step 04 - ENCODE hg38 blacklist v2
# Input:        ${BLACKLIST_URL}; one deduplicated BAM (chromosome naming)
# Output:       reference/hg38-blacklist.v2.bed
#               reference/hg38-blacklist.v2.bed3 (3 columns, BAM chromosome naming)
# Depends on:   curl, SAMtools
# Notes:        Needs network access; run after step 03:
#                 bash 04_prepare_blacklist.sh
# ============================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

if type module >/dev/null 2>&1; then
    module load samtools 2>/dev/null || module load SAMtools 2>/dev/null || true
fi
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 1; }

[ -s "${BLACKLIST_BED}" ] || curl -sSL "${BLACKLIST_URL}" | gunzip -c > "${BLACKLIST_BED}"
[ -s "${BLACKLIST_BED}" ] || { echo "ERROR: empty ${BLACKLIST_BED}" >&2; exit 1; }

BAM="$(ls "${DEDUP_DIR}"/*_dedup.bam 2>/dev/null | head -1)"
[ -n "${BAM}" ] || { echo "ERROR: no deduplicated BAM; run step 03" >&2; exit 1; }

match_chrom_to_bam "${BLACKLIST_BED}" "${BAM}" "${BLACKLIST_BED3}.tmp"
cut -f1-3 "${BLACKLIST_BED3}.tmp" | sort -k1,1 -k2,2n > "${BLACKLIST_BED3}"
rm -f "${BLACKLIST_BED3}.tmp"
echo "blacklist: $(wc -l < "${BLACKLIST_BED3}" | tr -d ' ') intervals -> ${BLACKLIST_BED3}"
