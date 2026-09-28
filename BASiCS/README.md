# Transcriptional noise quantification (BASiCS)

Quantification of transcriptional noise from single-cell RNA-sequencing data with
BASiCS. The analysis compares cell states in untreated cells (neuroblastoma and
hepatoblastoma models) and phases of a cisplatin treatment course.

For background on the method see the BASiCS repository
(https://github.com/catavallejos/BASiCS) and Vallejos et al. 2015,
*PLoS Comput Biol* (https://doi.org/10.1371/journal.pcbi.1004333).

## Script

`BASiCS_transcriptional_noise.R`, shown for MES cells across the cisplatin course.

1. **Data preparation.** Load the pre-processed scRNA-seq object
   (processing: https://github.com/PCEBrunaLab/Roux-et-al.-2024/tree/main/Single_Cell_RNA_Sequencing),
   add replicate information, remove unexpressed genes, and split cells into groups
   with an equal number of cells per group.
2. **Parameter estimation.** Estimate gene-level mean and over-dispersion with
   `BASiCS_MCMC` (N = 4000, Thin = 10, Burn = 2000, no spike-ins, replicates as batch).
3. **Differential over-dispersion.** Compare groups with `BASiCS_TestDE` (default parameters).
4. **Highly variable genes.** Detect highly variable genes within each group with
   `BASiCS_DetectHVG` (VarThreshold = 0.6; 0.3 for snSK-N-SH treated with DMSO, KDM5-C70 or T-5224).

## Usage

Set `INPUT_RDS` at the top of the script, then run:

```bash
Rscript BASiCS_transcriptional_noise.R
```

Chains are written to `chains/`; tables and plots to `results/`.

## Software

R (>= 4.2), BASiCS, SingleCellExperiment, Seurat, dplyr, ggplot2.
