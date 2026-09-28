# 1. snRNA-seq preprocessing
# Empty droplet removal, merging, normalisation and quality control.

library(SingleCellExperiment)
library(DropletUtils)
library(scuttle)
library(scater)
library(Seurat)
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(ggplot2)

setwd("/path/to/Multiome_project/")
set.seed(12)

## Load data ----
POT_data <- read10xCounts("cellranger_outputs/POT/raw_feature_bc_matrix/")
Cis_data <- read10xCounts("cellranger_outputs/Cisplatin_4weeksOFF/raw_feature_bc_matrix/")

## Remove empty droplets (FDR 0.1%) ----
#POT
calls <- emptyDrops(POT_data, assay.type = "counts", niters = 20000, ignore = 4999, lower = 500, retain = Inf)
POT.sce <- POT_data[, calls$FDR <= 0.001 & !is.na(calls$FDR)]
POT.sce$Sample <- "POT"

#Cisplatin 4 weeks OFF
calls <- emptyDrops(Cis_data, assay.type = "counts", niters = 20000, ignore = 4999, lower = 500, retain = Inf)
Cis.sce <- Cis_data[, calls$FDR <= 0.001 & !is.na(calls$FDR)]
Cis.sce$Sample <- "Cisplatin_4weeksOFF"

#Cell IDs as Sample#Barcode (same format as ArchR)
POT.sce$CellID <- paste0(POT.sce$Sample, "#", POT.sce$Barcode)
Cis.sce$CellID <- paste0(Cis.sce$Sample, "#", Cis.sce$Barcode)

## Merge samples and normalise ----
scRNA.sce <- cbind(POT.sce, Cis.sce)
colnames(scRNA.sce) <- scRNA.sce$CellID
scRNA.sce <- logNormCounts(scRNA.sce)

## Gene annotation ----
ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
ens.bm <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "external_gene_name", "description"),
                mart = ensembl, filters = "ensembl_gene_id", values = rownames(scRNA.sce))

#Mitochondrial, ribosomal and lincRNA genes
ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype == "lincRNA"]
linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]
mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]

#Add QC metrics
scRNA.sce <- addPerCellQCMetrics(scRNA.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))

#Keep named genes on main chromosomes, one row per gene name
ens.filt.bm <- ens.bm[ens.bm$chromosome_name %in% c(1:22, "X", "Y") & ens.bm$external_gene_name != "", ]
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
scRNA.sce <- scRNA.sce[ens.filt.bm$ensembl_gene_id, ]
rownames(scRNA.sce) <- ens.filt.bm$external_gene_name

## Cell cycle scoring ----
nb.seurat <- CreateSeuratObject(counts = counts(scRNA.sce), assay = "RNA", meta.data = as.data.frame(colData(scRNA.sce)))
nb.cycle.seurat <- NormalizeData(nb.seurat)
nb.cycle.seurat <- CellCycleScoring(nb.cycle.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
nb.seurat$Seurat.Phase <- nb.cycle.seurat$Phase
nb.seurat$Seurat.S <- nb.cycle.seurat$S.Score
nb.seurat$Seurat.G2M <- nb.cycle.seurat$G2M.Score

## Quality control filtering ----
#Plot QC metrics
p1 <- VlnPlot(nb.seurat, features = c("total", "subsets_mt_percent", "subsets_ribo_percent"), group.by = "Sample", ncol = 3, pt.size = 0)
ggsave(plot = p1, "plots/rna_qc_metrics.pdf", width = 12, height = 5)

#Genes: at least 3 reads, mitochondrial genes removed
rowcounts.filt <- rownames(nb.seurat)[Matrix::rowSums(GetAssayData(nb.seurat, layer = "counts")) >= 3]
genes.filt <- setdiff(rowcounts.filt, grep("^MT-", rownames(nb.seurat), value = TRUE))

#Cells: > 50 UMI, < 60% mitochondrial, > 0.25% ribosomal
nb.filt.seurat <- subset(nb.seurat,
                         cells = WhichCells(nb.seurat, expression = total > 50 & subsets_mt_percent < 60 & subsets_ribo_percent > 0.25),
                         features = genes.filt)
table(nb.filt.seurat$Sample)

## Save data files ----
saveRDS(nb.filt.seurat, "datafiles/nb_filt_seurat_qc.RDS")
