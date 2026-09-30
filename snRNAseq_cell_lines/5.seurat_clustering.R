# 5. Gene annotation, QC, cell cycle, SCTransform, PCA, UMAP and clustering
# Nuclei: >= 200 detected genes, < 25% mitochondrial, > 0.1% ribosomal reads.
# Genes: > 5 counts, mitochondrial genes removed.
# SCTransform (1,000 variable genes; cell cycle score and mitochondrial %
# regressed), PCA, UMAP and clustering on 30 PCs (resolution 0.2).
# Usage: Rscript 5.seurat_clustering.R <normalised.rds> <output.rds>

library(SingleCellExperiment)
library(scater)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]; out_file <- args[2]
set.seed(12345)
options(future.globals.maxSize = 12 * 1024^3)

nb.sce <- readRDS(in_file)

## Gene annotation: named genes on main chromosomes, one row per symbol ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ens.bm <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "external_gene_name"), mart = mart)
ens.filt <- ens.bm[ens.bm$chromosome_name %in% c(1:22, "X", "Y") & ens.bm$external_gene_name != "", ]
ens.filt <- ens.filt[!duplicated(ens.filt$external_gene_name), ]
idx <- match(rownames(nb.sce), ens.filt$ensembl_gene_id)
nb.sce <- nb.sce[!is.na(idx), ]
rownames(nb.sce) <- ens.filt$external_gene_name[idx[!is.na(idx)]]

## QC metrics ----
ens.86 <- genes(EnsDb.Hsapiens.v86)
linc.ids <- ens.86$gene_id[ens.86$gene_biotype == "lincRNA"]
linc.genes <- intersect(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% linc.ids], rownames(nb.sce))
ribo.genes <- grep("^RP[SL]", rownames(nb.sce), value = TRUE)
nb.sce <- addPerCellQCMetrics(nb.sce, flatten = TRUE, subsets = list(ribo = ribo.genes, linc = linc.genes))
nb.sce$subsets_mt_percent <- nb.sce$pct_mt      # mitochondrial % from step 3 (chrM is outside the kept chromosomes)
nb.sce$Condition <- factor(nb.sce$treatment, levels = c("DMSO", "C70", "T5"))

## Seurat object and cell cycle ----
meta <- as.data.frame(colData(nb.sce))[, c("Sample", "sample_id", "cell_line", "Condition", "replicate",
                                           "detected", "total", "subsets_mt_percent", "subsets_ribo_percent")]
nb.seurat <- CreateSeuratObject(counts = counts(nb.sce), assay = "RNA", meta.data = meta)
nb.seurat <- NormalizeData(nb.seurat, verbose = FALSE)
nb.seurat <- CellCycleScoring(nb.seurat, s.features = cc.genes.updated.2019$s.genes,
                              g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE)

## QC filtering ----
genes.keep <- setdiff(rownames(nb.sce)[Matrix::rowSums(counts(nb.sce)) > 5], grep("^MT-", rownames(nb.sce), value = TRUE))
cells.keep <- colnames(nb.sce)[nb.sce$detected >= 200 & nb.sce$subsets_mt_percent < 25 & nb.sce$subsets_ribo_percent > 0.1]
nb.seurat <- subset(nb.seurat, cells = cells.keep, features = genes.keep)
table(nb.seurat$Sample)

## SCTransform ----
nb.seurat$Seurat.Cycle.Score <- nb.seurat$S.Score - nb.seurat$G2M.Score
nb.seurat <- SCTransform(nb.seurat, method = "glmGamPoi", vst.flavor = "v2",
                         vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"),
                         variable.features.n = 1000, do.scale = TRUE, do.center = TRUE)

## PCA, UMAP, clustering ----
nb.seurat <- RunPCA(nb.seurat, npcs = 100, features = VariableFeatures(nb.seurat), verbose = FALSE)
nb.seurat <- RunUMAP(nb.seurat, reduction = "pca", dims = 1:30)
nb.seurat <- FindNeighbors(nb.seurat, reduction = "pca", dims = 1:30)
nb.seurat <- FindClusters(nb.seurat, resolution = 0.2)

DimPlot(nb.seurat, group.by = "seurat_clusters", label = TRUE, raster = TRUE)
ggsave("umap_clusters.pdf", width = 7, height = 6)
DimPlot(nb.seurat, group.by = "Condition", raster = TRUE)
ggsave("umap_condition.pdf", width = 7, height = 6)

saveRDS(nb.seurat, out_file)
