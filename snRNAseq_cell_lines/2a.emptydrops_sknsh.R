# 2a. Nucleus calling, SK-N-SH (emptyDrops)
# lower = 1000 UMI, niters = 10000, FDR <= 0.001.
# Usage: Rscript 2a.emptydrops_sknsh.R <raw_feature_bc_matrix.h5> <sample> <output.rds>

library(SingleCellExperiment)
library(DropletUtils)
library(Matrix)
library(BiocParallel)

args <- commandArgs(trailingOnly = TRUE)
h5_file <- args[1]; sample_id <- args[2]; out_file <- args[3]
register(SerialParam())

## Load raw matrix ----
sce <- read10xCounts(samples = h5_file, sample.names = sample_id, col.names = TRUE)
counts(sce) <- as(counts(sce), "dgCMatrix")
sce <- sce[, colSums(counts(sce)) > 0]

## Empty droplets ----
set.seed(42)
ed <- emptyDrops(counts(sce), lower = 1000, niters = 10000, retain = NULL,
                 test.ambient = FALSE, BPPARAM = SerialParam())
is_cell <- !is.na(ed$FDR) & ed$FDR <= 0.001
message(sample_id, ": ", sum(is_cell), " nuclei called")

sce <- sce[, is_cell]
sce$emptydrops_fdr <- ed$FDR[is_cell]
saveRDS(sce, out_file)
