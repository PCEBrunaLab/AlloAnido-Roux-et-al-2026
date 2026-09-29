# CUT&RUN processing analysis (SK-N-SH)

Scripts for the CUT&RUN section of the Methods: read processing, peak calling
on pooled replicates, CPM coverage tracks, and the comparison of AP-1-high and
AP-1-low promoters within the top H3K4me3 quartile.

- Targets: H3K4me3, KDM5A, ATF2, Fos, Jun
- Conditions (all at Awakening): Awakening (DMSO vehicle), C70 (KDM5-C70), T5 (T-5224)
- Two biological replicates per target and condition (30 libraries)

## Setup

1. Set `PROJ_DIR`, `GENOME_DIR` and `CUTTAG_DIR` in `00_config.sh`.
2. Place paired FASTQ files in `${PROJ_DIR}/raw_fastq/` as `<sample>_R1_001.fastq.gz` and `<sample>_R2_001.fastq.gz`.
3. Create `${PROJ_DIR}/samplesheet.csv` with the columns `sample_id,sample_name,condition,target,replicate`,
   and `${PROJ_DIR}/sample_list.txt` with one `sample_name` per line.
   Conditions are `Awakening`, `C70` and `T5`.
4. The parental SK-N-SH CUT&Tag H3K4me3 BAMs used in step 09 are expected at
   `${CUTTAG_DIR}/BAM_dedup/SKNSH_H3K4Me3_Rep{1,2}_dedup.bam`.

Keep all scripts in one directory. Each script sources `00_config.sh`.

## Run order

| Step | Script | Purpose |
|---|---|---|
| 01 | `01_trim_reads.sh` | Adapter and poly-G trimming (Cutadapt) |
| 02 | `02_map_reads.sh` | Alignment to GRCh38 (Bowtie2) |
| 03 | `03_remove_duplicates.sh` | Duplicate removal (Picard) |
| 04 | `04_prepare_blacklist.sh` | ENCODE hg38 blacklist v2 |
| 05 | `05_merge_replicates.sh` | Pooling of replicates (SAMtools) |
| 06 | `06_call_peaks.sh` | Peak calling on pooled BAMs (MACS3) |
| 07 | `07_peak_windows.sh` | Summit-centred windows per target |
| 08 | `08_cpm_tracks_heatmaps.sh` | CPM coverage tracks and occupancy heatmaps (deepTools) |
| 09 | `09_promoter_intervals.sh` | Promoter intervals from SK-N-SH CUT&Tag H3K4me3 |
| 10 | `10_promoter_ap1_split.sh` | AP-1-high vs AP-1-low promoters within the top H3K4me3 quartile |

```bash
sbatch 01_trim_reads.sh
sbatch 02_map_reads.sh
sbatch 03_remove_duplicates.sh
bash   04_prepare_blacklist.sh        # requires network access
sbatch 05_merge_replicates.sh
sbatch 06_call_peaks.sh
sbatch 07_peak_windows.sh
sbatch 08_cpm_tracks_heatmaps.sh
sbatch 09_promoter_intervals.sh
sbatch 10_promoter_ap1_split.sh
```

## Parameters

| Tool | Version | Settings |
|---|---|---|
| Cutadapt | 4.9 | TruSeq adapters; `-O 1 --nextseq-trim=20 -q 0 -m 20` |
| Bowtie2 | 2.4.2 | `--local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700` |
| Picard | 2.23.8 | `MarkDuplicates REMOVE_DUPLICATES=true`; no mapping-quality filter |
| SAMtools | 1.11 | deduplicated replicates merged per condition and target |
| MACS3 | 3.0.4 | `callpeak -f BAMPE -g hs -q 0.1 --nomodel --keep-dup all`; no control; local lambda |
| deepTools | 3.5.5 | `bamCoverage --normalizeUsing CPM --binSize 50 --extendReads`, blacklist excluded; `computeMatrix reference-point` ±3 kb, 50 bp bins |
| bedtools | 2.29.2 | blacklist and interval operations |

**Promoter intervals (step 09).** The two SK-N-SH CUT&Tag H3K4me3 replicates are merged and called
with the same MACS3 settings. Peaks on the main chromosomes that do not overlap the blacklist
form the promoter set.

**AP-1 split (step 10).**
1. Mean CPM of H3K4me3, KDM5A, ATF2, Fos and Jun is computed over each promoter interval (`multiBigwigSummary`).
2. Ranking scores are averaged over the three conditions (Awakening, C70, T5).
3. The top 25% of promoters by H3K4me3 are kept.
4. Within these, the top and bottom 25% by ATF2, Fos or Jun define the AP-1-high and AP-1-low groups.
5. KDM5A and H3K4me3 are compared in each condition as the ratio of median CPM (high / low).
   The H3K4me3 ratio is reported as a balance check.

Heatmap colour scales are set to the 98th percentile of each panel.

## Outputs

| Path | Content |
|---|---|
| `results/peaks/<Condition>/<Target>/` | MACS3 peaks and summits |
| `results/bigwig/<Condition>_<Target>.cpm.bw` | CPM coverage tracks |
| `results/heatmaps/` | Occupancy heatmaps per target |
| `results/promoters/promoter.bed` | Promoter intervals |
| `results/promoter_ap1_split/results.tsv` | KDM5A and H3K4me3 ratios per factor |
| `results/promoter_ap1_split/<Factor>/` | `summary.tsv`, `regions/`, heatmaps |

## Software

Cutadapt 4.9, Bowtie2 2.4.2, Picard 2.23.8, SAMtools 1.11, MACS3 3.0.4, bedtools 2.29.2,
deepTools 3.5.5, Python 3 with numpy. Scripts are written for SLURM; `#SBATCH` resources may need
adjusting for other clusters.

## Software environment

Conda environments with the tool versions used are in `envs/`. Create one with:

```bash
mamba env create -f envs/<file>.yml
```

| File | Used for |
|---|---|
| `envs/cutrun_processing.yml` | Steps 01-05 |
| `envs/cutrun_analysis.yml` | Steps 06-10 |
