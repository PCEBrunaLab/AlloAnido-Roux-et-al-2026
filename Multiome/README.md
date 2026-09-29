# Single-nuclei multiomics (snRNA-seq and snATAC-seq)

Scripts for the single-nuclei multiome analysis of SK-N-SH cells before cisplatin
treatment (`POT`) and after recovery (`Cisplatin_4weeksOFF`).

Raw FASTQs were processed with 10x Genomics `cellranger-arc count` (v2.0.2, GRCh38-2024-A).
Each script expects the cellranger-arc outputs in `cellranger_outputs/<sample>/` and writes to
`datafiles/` and `plots/` in the project folder set by `setwd()` at the top of the script.

## Scripts

| Script | Methods section |
|---|---|
| `1.RNA_Preprocessing.R` | Single nuclei RNA preprocessing: empty droplets, merging, normalisation, QC |
| `2.RNA_Normalisation_Clustering.R` | SCTransform, PCA, UMAP, clustering, cluster markers, gene set enrichment |
| `3.RNA_Cell_State_Scores.R` | Cell state scores (ADRN / MES) |
| `4.ATAC_Preprocessing.R` | Single nuclei ATAC preprocessing: Arrow files, QC, doublets, LSI, UMAP |
| `5.ATAC_Peak_Calling.R` | Single nuclei ATAC peak calling (Signac) and combined peak set |
| `6.ATAC_Differential_Peaks_Motifs.R` | Differential accessibility and transcription factor motif enrichment |

## Software

R with DropletUtils, scuttle, scater, Seurat, glmGamPoi, clusterProfiler, msigdbr,
Signac, ArchR (v1.0.3), GenomicRanges, chromVARmotifs and MACS.

## Software environment

Conda environments with the tool versions used are in `envs/`. Create one with:

```bash
mamba env create -f envs/<file>.yml
```

| File | Used for |
|---|---|
| `envs/multiome_R.yml` | All six scripts (ArchR installed from GitHub, see the file header) |
| `envs/multiome_macs.yml` | MACS3 for peak calling in script 5 |
