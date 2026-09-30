# 3. Nucleus QC and doublet removal (per sample)
# 500-10,000 detected genes, >= 1000 UMI, <= 25% mitochondrial reads,
# then scDblFinder (expected doublet rate = 0.008 x nuclei / 1000).
# Input: emptyDrops nuclei (SK-N-SH, .rds) or CellBender filtered H5 (HuH6, .h5)
# Usage: Rscript 3.filter_nuclei_qc.R <input.rds|input.h5> <sample> <output.rds>

library(SingleCellExperiment)
library(DropletUtils)
library(scDblFinder)
library(Matrix)
library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
in_file <- args[1]; sample_id <- args[2]; out_file <- args[3]
register(SerialParam())
set.seed(42)

## Load nuclei ----
if (grepl("\\.h5$", in_file)) {
  sce <- read10xCounts(in_file, col.names = TRUE)      # CellBender filtered matrix
  sce$Sample <- sample_id
} else {
  sce <- readRDS(in_file)                              # emptyDrops nuclei
}
counts(sce) <- as(counts(sce), "dgCMatrix")

## Per-cell metrics ----
symbols <- as.character(rowData(sce)$Symbol)
is_mt <- grepl("^MT-", symbols)
sce$n_genes <- colSums(counts(sce) > 0)
sce$n_umi   <- colSums(counts(sce))
sce$pct_mt  <- 100 * colSums(counts(sce)[is_mt, , drop = FALSE]) / sce$n_umi

## Filters ----
keep <- sce$n_genes >= 500 & sce$n_genes <= 10000 & sce$n_umi >= 1000 & sce$pct_mt <= 25
sce <- sce[, keep]

## Doublets ----
set.seed(42)
sce <- scDblFinder(sce, dbr = 0.008 * ncol(sce) / 1000, BPPARAM = SerialParam())
sce <- sce[, sce$scDblFinder.class == "singlet"]
message(sample_id, ": ", ncol(sce), " nuclei after QC")

saveRDS(sce, out_file)
