# 4. Normalisation across samples (scran / batchelor), per cell line
# Deconvolution size factors within each sample, then multiBatchNorm across
# samples. Genes for size factors: detected in >= 3 nuclei and mean > 0.01.
# Usage: Rscript 4.normalise.R <SK-N-SH|HuH6> <qc_dir> <output.rds>

library(SingleCellExperiment)
library(scran)
library(batchelor)
library(Matrix)
library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
cell_line <- args[1]; qc_dir <- args[2]; out_file <- args[3]
set.seed(42)
bp <- MulticoreParam(workers = 4)

sheet <- read.delim("sample_sheet.tsv", stringsAsFactors = FALSE)
sheet <- sheet[sheet$cell_line == cell_line, ]

## Load QC'd samples and add metadata ----
sce.list <- lapply(sheet$condition_rep, function(s) {
  x <- readRDS(file.path(qc_dir, paste0(s, "_SCE_qc.RDS")))
  rowData(x) <- rowData(x)[, !grepl("^scDblFinder", colnames(rowData(x))), drop = FALSE]
  colnames(x) <- paste(s, x$Barcode, sep = "_")
  meta <- sheet[sheet$condition_rep == s, ]
  x$Sample <- s; x$sample_id <- meta$sample_id; x$cell_line <- meta$cell_line
  x$treatment <- meta$treatment; x$replicate <- as.character(meta$replicate)
  x
})
names(sce.list) <- sheet$condition_rep

## Genes used for size factors ----
all.counts <- do.call(cbind, lapply(sce.list, counts))
keep_sf <- rowSums(all.counts > 0) >= 3 & rowMeans(all.counts) > 0.01
rm(all.counts)

## Size factors within each sample ----
for (s in names(sce.list)) {
  cl <- quickCluster(sce.list[[s]], subset.row = keep_sf, method = "igraph", min.size = 100, BPPARAM = bp)
  sce.list[[s]] <- computeSumFactors(sce.list[[s]], clusters = cl, subset.row = keep_sf,
                                     positive = TRUE, BPPARAM = bp)
}

## Scale across samples ----
sce <- do.call(cbind, do.call(multiBatchNorm, c(unname(sce.list), list(min.mean = 0.01, normalize.all = TRUE))))
saveRDS(sce, out_file)
