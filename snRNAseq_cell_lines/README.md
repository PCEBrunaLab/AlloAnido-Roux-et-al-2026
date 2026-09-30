# snRNA-seq of SK-N-SH and HuH6 under KDM5 and AP-1 inhibition

Single-nucleus RNA-seq of SK-N-SH (neuroblastoma) and HuH6 (hepatoblastoma) cells treated
with DMSO, KDM5-C70 or T-5224, two replicates each (12 libraries, 10x Genomics Single Cell 3' v4).

## Scripts

Run in order. `sample_sheet.tsv` maps each library to cell line, treatment and replicate.

| Script | Step |
|---|---|
| `1.cellranger_count.sh` | Cell Ranger 8.0.0 `count`, reference `refdata-gex-GRCh38-2020-A` |
| `2a.emptydrops_sknsh.R` | Nucleus calling, SK-N-SH: emptyDrops (lower = 1,000 UMI, FDR ≤ 0.001) |
| `2b.cellbender_huh6.sh` | Nucleus calling and ambient RNA removal, HuH6: CellBender 0.4.0 (FPR 0.01, 150 epochs) |
| `3.filter_nuclei_qc.R` | 500–10,000 genes, ≥ 1,000 UMI, ≤ 25% mitochondrial; scDblFinder doublet removal |
| `4.normalise.R` | Deconvolution size factors per sample, multiBatchNorm across samples (per cell line) |
| `5.seurat_clustering.R` | ≥ 200 genes, < 25% mitochondrial, > 0.1% ribosomal; SCTransform (1,000 genes, cell cycle and mitochondrial % regressed); PCA, UMAP and clustering (30 PCs, resolution 0.2) |
| `6.cell_states_sknsh.R` | ADRN / MES scores and per-nucleus cell state (van Groningen et al. 2017), SK-N-SH |

```bash
sbatch 1.cellranger_count.sh
Rscript 2a.emptydrops_sknsh.R cellranger/SKNSH_DMSO_1/outs/raw_feature_bc_matrix.h5 SKNSH_DMSO_1 SKNSH_DMSO_1_SCE.RDS   # per SK-N-SH sample
sbatch 2b.cellbender_huh6.sh
Rscript 3.filter_nuclei_qc.R SKNSH_DMSO_1_SCE.RDS SKNSH_DMSO_1 qc/SKNSH_DMSO_1_SCE_qc.RDS                              # per sample
Rscript 3.filter_nuclei_qc.R cellbender/HuH6_DMSO_1/HuH6_DMSO_1_cellbender_filtered.h5 HuH6_DMSO_1 qc/HuH6_DMSO_1_SCE_qc.RDS
Rscript 4.normalise.R SK-N-SH qc sknsh_normalised.rds
Rscript 5.seurat_clustering.R sknsh_normalised.rds sknsh_seurat.rds
Rscript 6.cell_states_sknsh.R sknsh_seurat.rds vanGroningen_2017.xlsx sknsh_seurat_states.rds
```

Repeat steps 4 and 5 with `HuH6` for the HuH6 samples.

## Software environment

| File | Used for |
|---|---|
| `envs/snrnaseq_R.yml` | Scripts 2a and 3–6 |
| `envs/cellbender.yml` | Script 2b |

Cell Ranger 8.0.0 is installed separately from 10x Genomics.
