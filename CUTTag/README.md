# CUT&Tag Analysis Pipeline

A complete bioinformatics pipeline for processing and analyzing CUT&Tag (Cleavage Under Targets and Tagmentation) data, from raw sequencing reads to differential binding analysis.

## Overview

This pipeline was developed for analyzing CUT&Tag data profiling histone modifications (H3K4me3) and chromatin regulators (KDM5A, KDM5B) across treatment time points in neuroblastoma cells. It can be adapted for any CUT&Tag experiment with paired-end sequencing data.

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 00 | `00_config.sh` | Central configuration (paths, parameters, sample mapping) |
| 01 | `01_mapping.sh` | Bowtie2 alignment with CUT&Tag-specific parameters |
| 02 | `02_quality_control.sh` | QC metrics: alignment rate, duplication, library complexity |
| 03 | `03_remove_duplicates.sh` | PCR duplicate removal (Picard MarkDuplicates) |
| 04a | `04a_peak_calling.sh` | MACS2 peak calling (treatment vs. untreated control) |
| 04b | `04b_peak_calling_untreated.sh` | MACS2 peak calling (untreated vs. IgG control) |
| 05 | `05_frip_analysis.R` | FRiP, peak numbers, peak width analysis |
| 06 | `06_bigwig_generation.sh` | RPGC-normalized BigWig files (deepTools bamCoverage) |
| 07 | `07_heatmaps.sh` | Genome-wide TSS and peak-centered heatmaps (deepTools) |
| 08 | `08_gene_specific_analysis.sh` | Signal analysis at custom gene lists |
| 09 | `09_cellstate_analysis.sh` | Cell state signature analysis (ADRN/MES/Intermediate) |
| 10 | `10_peak_annotation.R` | Genomic distribution analysis (ChIPseeker) |
| 11 | `11_differential_binding.R` | Differential binding & pathway enrichment (5 databases) |

## Quick Start

1. Clone this repository
2. Edit `00_config.sh` with your paths and sample mapping
3. Run scripts sequentially (01 → 11)

```bash
# Edit configuration
vim 00_config.sh

# Run pipeline
sbatch 01_mapping.sh
sbatch 02_quality_control.sh
sbatch 03_remove_duplicates.sh
sbatch 04a_peak_calling.sh
sbatch 04b_peak_calling_untreated.sh
Rscript 05_frip_analysis.R
sbatch 06_bigwig_generation.sh
sbatch 07_heatmaps.sh
sbatch 08_gene_specific_analysis.sh
sbatch 09_cellstate_analysis.sh
Rscript 10_peak_annotation.R
Rscript 11_differential_binding.R
```

## Dependencies

### Command-line tools
- Bowtie2 (≥ 2.4)
- SAMtools (≥ 1.11)
- Picard Tools
- MACS2 (≥ 2.2.7)
- deepTools (≥ 3.5) — bamCoverage, computeMatrix, plotHeatmap, plotProfile
- BEDTools (≥ 2.29)
- HOMER (mergePeaks)

### R packages
- tidyverse, ggplot2, ggpubr
- GenomicRanges, Rsamtools
- ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db
- clusterProfiler, ReactomePA, enrichplot, msigdbr
- pheatmap, VennDiagram, ComplexHeatmap

## Experimental Design

The pipeline supports multi-condition, multi-target CUT&Tag experiments with biological replicates:

- **Conditions**: Treatment time points (e.g., Untreated → Day 7 → Day 14 → Recovery)
- **Targets**: Histone marks and chromatin regulators (e.g., H3K4me3, KDM5A, KDM5B)
- **Controls**: IgG negative control, untreated time-point control
- **Replicates**: Biological replicates processed separately and merged

## Output Structure

```
project/
├── bam_files/              # Mapped BAMs (Step 01)
├── bam_files_dedup/        # Deduplicated BAMs (Step 03)
├── QC_final/               # QC summary and plots (Step 02)
├── peakCalling/MACS2/      # Peak files (Steps 04a, 04b)
├── alignment/bigwig/       # Normalized BigWig files (Step 06)
├── heatmap/                # Heatmap outputs (Steps 07-09)
│   ├── final_timeseries/   # Genome-wide heatmaps
│   ├── gene_specific_peaks/# Gene-list specific analysis
│   └── cellstate_genelist/ # Cell state analysis
└── Analysis/               # Annotation and enrichment (Steps 10-11)
    ├── Genomic_Distribution_Detailed/
    └── Differential_Binding_Pathways/
```

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| MACS2 q-value | 0.1 | Peak calling significance threshold |
| Peak format | BAMPE | Paired-end BAM format for MACS2 |
| BigWig normalization | RPGC | Reads Per Genomic Content |
| Heatmap distance | 3000 bp | Upstream/downstream of reference point |
| Peak-gene distance | 5000 bp | Maximum distance for peak-to-gene assignment |

## Citation

If you use this pipeline, please cite the relevant tools:
- MACS2: Zhang et al., Genome Biology 2008
- deepTools: Ramirez et al., Nucleic Acids Research 2016
- ChIPseeker: Yu et al., Bioinformatics 2015
- clusterProfiler: Wu et al., The Innovation 2021
