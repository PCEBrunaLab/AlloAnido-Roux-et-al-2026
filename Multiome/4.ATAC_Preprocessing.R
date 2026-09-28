# 4. snATAC-seq preprocessing and visualisation (ArchR)
# Arrow files, QC, doublet removal, iterative LSI, clustering and UMAP.

library(ArchR)

setwd("/path/to/Multiome_project/")
addArchRThreads(threads = 4)
addArchRGenome("hg38")
set.seed(12)

## Create Arrow files ----
#Cells kept with > 1000 fragments and TSS enrichment > 4
ArrowFiles <- createArrowFiles(inputFiles = c("cellranger_outputs/POT/atac_fragments.tsv.gz",
                                              "cellranger_outputs/Cisplatin_4weeksOFF/atac_fragments.tsv.gz"),
                               sampleNames = c("POT", "Cisplatin_4weeksOFF"),
                               minTSS = 4,
                               minFrags = 1000,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               force = TRUE)

## Infer doublets ----
doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)

## Create ArchR project ----
project <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "ArchR_project/", copyArrows = TRUE)

#Plot number of fragments vs TSS enrichment
df <- getCellColData(project, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(x = df[, 1], y = df[, 2], colorDensity = TRUE, continuousSet = "sambaNight",
             xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment") +
  geom_hline(yintercept = 4, lty = "dashed") +
  geom_vline(xintercept = 3, lty = "dashed")
ggsave(plot = p, "plots/atac_QC_metrics.pdf", width = 7, height = 7)

## Filter doublets ----
project <- filterDoublets(project, filterRatio = 1.2)

## Dimensionality reduction (iterative LSI, TF-IDF) ----
project <- addIterativeLSI(ArchRProj = project,
                           useMatrix = "TileMatrix",
                           name = "LSI_ATAC",
                           iterations = 2,
                           clusterParams = list(resolution = c(0.2), n.start = 10),
                           varFeatures = 25000,
                           dimsToUse = 1:30,
                           force = TRUE)

## Clustering and UMAP ----
project <- addClusters(input = project, reducedDims = "LSI_ATAC", method = "Seurat",
                       name = "Clusters_LSI_ATAC", resolution = 0.6, force = TRUE)

project <- addUMAP(ArchRProj = project, reducedDims = "LSI_ATAC", name = "UMAP_LSI",
                   nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)

pdf("plots/atac_umap.pdf")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP_LSI")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters_LSI_ATAC", embedding = "UMAP_LSI")
dev.off()

## Save ----
saveRDS(getCellColData(project), "datafiles/atac_metadata.RDS")
saveRDS(getEmbedding(ArchRProj = project, embedding = "UMAP_LSI", returnDF = TRUE), "datafiles/UMAP_LSI_coordinates.RDS")
saveArchRProject(ArchRProj = project, outputDirectory = "ArchR_project/", load = FALSE)
