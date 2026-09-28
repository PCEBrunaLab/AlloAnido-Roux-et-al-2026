#!/usr/bin/env Rscript
#==================================================================================
# Step 11: Differential Binding & Pathway Analysis
# Identifies condition-specific peaks and runs enrichment analysis
#
# Input:  Peak files from Step 04 (narrowPeak)
# Output: Differential peak sets, gene lists, enrichment results (5 databases)
#
# Databases: GO-BP, KEGG, Reactome, MSigDB Hallmark, MSigDB C6 Oncogenic
# Dependencies: ChIPseeker, clusterProfiler, ReactomePA, msigdbr, VennDiagram
#==================================================================================

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(enrichplot)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(VennDiagram)
  library(RColorBrewer)
  library(pheatmap)
  library(tibble)
  library(rtracklayer)
})

# Setup
projPath <- Sys.getenv("PROJ_DIR", "/path/to/your/project")
peakDir  <- file.path(projPath, "peakCalling/MACS2/human/Untreated_control")
outDir   <- file.path(projPath, "Analysis/Differential_Binding_Pathways")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
conditions <- c("D7_Cis", "D14_Cis", "D14_C70")
replicates <- c("rep1", "rep2")

# Prepare MSigDB gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

c6_sets <- msigdbr(species = "Homo sapiens", category = "C6")
c6_list <- split(c6_sets$gene_symbol, c6_sets$gs_name)

# Load KDM5A peaks (primary target for differential analysis)
# NOTE: Change "KDM5A" to your marker of interest
TARGET_MARKER <- "KDM5A"

cat(sprintf("Loading %s peaks...\n", TARGET_MARKER))
marker_peaks <- list()

for(condition in conditions) {
  all_peaks <- GRanges()
  for(rep in replicates) {
    peak_file <- file.path(peakDir, condition, TARGET_MARKER, rep,
                           "macs2_peak_q0.1_peaks.narrowPeak.bed")
    if(file.exists(peak_file)) {
      peaks <- read.table(peak_file, sep="\t", stringsAsFactors = FALSE)
      if(ncol(peaks) >= 3 && nrow(peaks) > 0) {
        colnames(peaks)[1:3] <- c("chr", "start", "end")
        if(!grepl("^chr", peaks$chr[1])) peaks$chr <- paste0("chr", peaks$chr)
        peaks <- peaks[peaks$start >= 0 & peaks$end > peaks$start, ]
        if(nrow(peaks) > 0) {
          all_peaks <- c(all_peaks, GRanges(seqnames = peaks$chr,
                        ranges = IRanges(start = peaks$start, end = peaks$end)))
        }
      }
    }
  }
  if(length(all_peaks) > 0) {
    marker_peaks[[condition]] <- reduce(all_peaks)
    cat(sprintf("  %s: %d peaks\n", condition, length(marker_peaks[[condition]])))
  }
}

# Identify differential peak sets
cat("\n=== DIFFERENTIAL PEAK SETS ===\n")

all_D7 <- marker_peaks[["D7_Cis"]]
all_D14 <- marker_peaks[["D14_Cis"]]
all_C70 <- marker_peaks[["D14_C70"]]

# D14-specific (awakening)
overlaps_D14_D7 <- findOverlaps(all_D14, all_D7)
D14_specific <- all_D14[-queryHits(overlaps_D14_D7)]
D7_D14_shared <- all_D14[queryHits(overlaps_D14_D7)]

# C70 comparison
overlaps_D14_C70 <- findOverlaps(all_D14, all_C70)
C70_lost <- all_D14[-queryHits(overlaps_D14_C70)]
C70_maintained <- all_D14[queryHits(overlaps_D14_C70)]

peak_sets <- list(
  "D7_All" = all_D7, "D14_Specific" = D14_specific,
  "D7_D14_Shared" = D7_D14_shared, "C70_Lost" = C70_lost,
  "C70_Maintained" = C70_maintained, "D14_All" = all_D14, "C70_All" = all_C70
)

for(name in names(peak_sets)) {
  cat(sprintf("  %s: %d peaks\n", name, length(peak_sets[[name]])))
  export(peak_sets[[name]], file.path(outDir, paste0(TARGET_MARKER, "_", name, "_peaks.bed")), format = "bed")
}

# Annotate peaks to genes
cat("\n=== ANNOTATING PEAKS ===\n")
gene_lists <- list()

for(name in names(peak_sets)) {
  if(length(peak_sets[[name]]) == 0) next
  anno <- annotatePeak(peak_sets[[name]], tssRegion = c(-3000, 1000),
                       TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  genes <- unique(as.data.frame(anno)$SYMBOL[!is.na(as.data.frame(anno)$SYMBOL)])
  gene_lists[[name]] <- genes
  write.table(genes, file.path(outDir, paste0("Genes_", name, ".txt")),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Enrichment analysis function
run_enrichment <- function(genes, name) {
  if(length(genes) < 10) return(NULL)
  gene_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if(nrow(gene_ids) < 10) return(NULL)
  entrez_ids <- gene_ids$ENTREZID
  results <- list()

  tryCatch({
    go <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
    if(!is.null(go) && nrow(go) > 0) results$GO_BP <- go
  }, error = function(e) NULL)

  tryCatch({
    kegg <- enrichKEGG(entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)
    if(!is.null(kegg) && nrow(kegg) > 0) {
      results$KEGG <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
  }, error = function(e) NULL)

  tryCatch({
    react <- enrichPathway(entrez_ids, pvalueCutoff = 0.05, readable = TRUE, organism = "human")
    if(!is.null(react) && nrow(react) > 0) results$Reactome <- react
  }, error = function(e) NULL)

  tryCatch({
    h <- enricher(genes, TERM2GENE = data.frame(
      term = rep(names(hallmark_list), sapply(hallmark_list, length)),
      gene = unlist(hallmark_list)), pvalueCutoff = 0.05)
    if(!is.null(h) && nrow(h) > 0) results$Hallmark <- h
  }, error = function(e) NULL)

  tryCatch({
    c6 <- enricher(genes, TERM2GENE = data.frame(
      term = rep(names(c6_list), sapply(c6_list, length)),
      gene = unlist(c6_list)), pvalueCutoff = 0.05)
    if(!is.null(c6) && nrow(c6) > 0) results$C6_Oncogenic <- c6
  }, error = function(e) NULL)

  return(results)
}

# Run enrichment
cat("\n=== PATHWAY ENRICHMENT ===\n")
enrichment_results <- list()
for(set_name in names(gene_lists)) {
  enrichment_results[[set_name]] <- run_enrichment(gene_lists[[set_name]], set_name)
}

# Save results and create plots
for(set_name in names(enrichment_results)) {
  if(is.null(enrichment_results[[set_name]])) next
  set_dir <- file.path(outDir, "Enrichment", set_name)
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)
  for(db in names(enrichment_results[[set_name]])) {
    result <- enrichment_results[[set_name]][[db]]
    if(!is.null(result) && nrow(result) > 0) {
      write.csv(as.data.frame(result), file.path(set_dir, paste0(db, ".csv")), row.names = FALSE)
      tryCatch({
        ggsave(file.path(set_dir, paste0(db, "_dotplot.pdf")),
               dotplot(result, showCategory = 20), width = 12, height = 10)
      }, error = function(e) NULL)
    }
  }
}

# Venn diagram
pdf(file.path(outDir, "Gene_Overlap_Venn.pdf"), width = 8, height = 8)
venn.plot <- venn.diagram(
  x = list("D7" = gene_lists[["D7_All"]],
           "D14" = gene_lists[["D14_All"]],
           "C70" = gene_lists[["C70_All"]]),
  filename = NULL,
  fill = c("#E41A1C", "#4DAF4A", "#984EA3"), alpha = 0.5,
  main = paste(TARGET_MARKER, "Target Gene Overlap"))
grid.draw(venn.plot)
dev.off()

cat("\nAnalysis complete! Results in:", outDir, "\n")
