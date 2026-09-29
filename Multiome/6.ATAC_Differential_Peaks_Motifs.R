# 6. Chromatin accessibility and transcription factor motif analysis
# Cells with both snATAC-seq and snRNA-seq; differential peaks between
# pre- and post-treatment cells of the same state; CIS-BP motif enrichment.

library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("/path/to/Multiome_project/")
addArchRThreads(threads = 4)
addArchRGenome("hg38")
set.seed(12)

project <- loadArchRProject("ArchR_project/")
nb.seurat <- readRDS("datafiles/nb_seurat_AMT.rds")

## Keep cells with both ATAC and RNA ----
common_cells <- intersect(getCellNames(project), colnames(nb.seurat))
project <- project[common_cells, ]

#Add RNA cell state to ArchR project
project$AMT.state <- as.character(nb.seurat$AMT.state[getCellNames(project)])
project$Condition_AMT <- paste0(project$Sample, "_", project$AMT.state)
table(project$Condition_AMT)

#Overlay RNA cell state on the ATAC UMAP
pdf("plots/RNA_cell_state_on_ATAC_UMAP.pdf")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "AMT.state", embedding = "UMAP_LSI")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP_LSI")
dev.off()

## Differential peaks: pre vs post treatment ----
#Positive Log2FC = higher accessibility in POT (useGroups)
diffPeaks_MES <- getMarkerFeatures(ArchRProj = project,
                                   useMatrix = "PeakMatrix",
                                   groupBy = "Condition_AMT",
                                   bias = c("TSSEnrichment", "log10(nFrags)"),
                                   testMethod = "wilcoxon",
                                   normBy = "ReadsInTSS",
                                   useGroups = "POT_MES",
                                   bgdGroups = "Cisplatin_4weeksOFF_MES",
                                   maxCells = 2100)

diffPeaks_ADRN <- getMarkerFeatures(ArchRProj = project,
                                    useMatrix = "PeakMatrix",
                                    groupBy = "Condition_AMT",
                                    bias = c("TSSEnrichment", "log10(nFrags)"),
                                    testMethod = "wilcoxon",
                                    normBy = "ReadsInTSS",
                                    useGroups = "POT_ADRN",
                                    bgdGroups = "Cisplatin_4weeksOFF_ADRN",
                                    maxCells = 3500)

write.csv(as.data.frame(getMarkers(diffPeaks_MES, cutOff = "FDR <= 1")[[1]]), "datafiles/diffPeaks_MES.csv")
write.csv(as.data.frame(getMarkers(diffPeaks_ADRN, cutOff = "FDR <= 1")[[1]]), "datafiles/diffPeaks_ADRN.csv")

#Volcano plots (FDR <= 0.1, abs(Log2FC) >= 1)
pdf("plots/Volcano_plots_differential_peaks.pdf")
plotMarkers(seMarker = diffPeaks_MES, name = "POT_MES", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotMarkers(seMarker = diffPeaks_ADRN, name = "POT_ADRN", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
dev.off()

## Motif enrichment (CIS-BP, chromVARmotifs human_pwms_v2) ----
project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif", force = TRUE)

#Rank motifs by enrichment, one row per TF
rank_motifs <- function(motifs){
  df <- data.frame(TF = rownames(motifs), enrichment = assay(motifs, "Enrichment")[, 1])
  df <- df[order(df$enrichment, decreasing = TRUE), ]
  df$TF <- gsub("_[0-9]+$", "", df$TF)
  df <- df[!duplicated(df$TF), ]
  df$rank <- seq_len(nrow(df))
  df
}

#MES
motifs_MES_POT <- peakAnnoEnrichment(seMarker = diffPeaks_MES, ArchRProj = project, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
motifs_MES_Cis <- peakAnnoEnrichment(seMarker = diffPeaks_MES, ArchRProj = project, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
write.csv(rank_motifs(motifs_MES_POT), "datafiles/MES_motifs_higher_in_POT.csv", row.names = FALSE)
write.csv(rank_motifs(motifs_MES_Cis), "datafiles/MES_motifs_higher_in_Cisplatin_4weeksOFF.csv", row.names = FALSE)

#ADRN
motifs_ADRN_POT <- peakAnnoEnrichment(seMarker = diffPeaks_ADRN, ArchRProj = project, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
motifs_ADRN_Cis <- peakAnnoEnrichment(seMarker = diffPeaks_ADRN, ArchRProj = project, peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
write.csv(rank_motifs(motifs_ADRN_POT), "datafiles/ADRN_motifs_higher_in_POT.csv", row.names = FALSE)
write.csv(rank_motifs(motifs_ADRN_Cis), "datafiles/ADRN_motifs_higher_in_Cisplatin_4weeksOFF.csv", row.names = FALSE)

## Save ----
saveArchRProject(ArchRProj = project, outputDirectory = "ArchR_project/", load = FALSE)
