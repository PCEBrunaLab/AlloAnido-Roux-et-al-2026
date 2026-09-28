# 5. snATAC-seq peak calling (Signac)
# Peaks called per sample, blacklist and sex chromosome peaks removed,
# then combined with GenomicRanges::reduce and added to the ArchR project.

library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ArchR)

setwd("/path/to/Multiome_project/")
set.seed(12)
macs_path <- "/path/to/macs"

#Gene annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

#Cells passing ArchR QC (Sample#Barcode)
archr_cells <- rownames(readRDS("datafiles/atac_metadata.RDS"))

## Peak calling POT ----
counts <- Read10X_h5("cellranger_outputs/POT/filtered_feature_bc_matrix.h5")[["Peaks"]]
cells <- intersect(colnames(counts), gsub("POT#", "", grep("^POT#", archr_cells, value = TRUE)))
sample_atac <- CreateSeuratObject(counts = CreateChromatinAssay(counts = counts[, cells], sep = c(":", "-"),
                                                                fragments = "cellranger_outputs/POT/atac_fragments.tsv.gz",
                                                                annotation = annotation),
                                  assay = "ATAC")
POT.peaks <- CallPeaks(object = sample_atac, macs2.path = macs_path)

## Peak calling Cisplatin 4 weeks OFF ----
counts <- Read10X_h5("cellranger_outputs/Cisplatin_4weeksOFF/filtered_feature_bc_matrix.h5")[["Peaks"]]
cells <- intersect(colnames(counts), gsub("Cisplatin_4weeksOFF#", "", grep("^Cisplatin_4weeksOFF#", archr_cells, value = TRUE)))
sample_atac <- CreateSeuratObject(counts = CreateChromatinAssay(counts = counts[, cells], sep = c(":", "-"),
                                                                fragments = "cellranger_outputs/Cisplatin_4weeksOFF/atac_fragments.tsv.gz",
                                                                annotation = annotation),
                                  assay = "ATAC")
Cis.peaks <- CallPeaks(object = sample_atac, macs2.path = macs_path)

## Remove blacklist, non-standard and sex chromosome peaks ----
filter_peaks <- function(peaks){
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  peaks[seqnames(peaks) %in% paste0("chr", 1:22)]
}
POT.peaks <- filter_peaks(POT.peaks)
Cis.peaks <- filter_peaks(Cis.peaks)

## Combine peaks across samples ----
combined.peaks <- reduce(x = c(POT.peaks, Cis.peaks))
combined.peaks
saveRDS(combined.peaks, "datafiles/combined_peaks_reduced.RDS")

## Add peak set to ArchR project ----
addArchRThreads(threads = 1)
addArchRGenome("hg38")
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE", RHDF5_USE_FILE_LOCKING = "FALSE")

project <- loadArchRProject("ArchR_project/")
project <- addPeakSet(ArchRProj = project, peakSet = combined.peaks, force = TRUE)
project <- addPeakMatrix(project)
saveArchRProject(ArchRProj = project, outputDirectory = "ArchR_project/", load = FALSE)
