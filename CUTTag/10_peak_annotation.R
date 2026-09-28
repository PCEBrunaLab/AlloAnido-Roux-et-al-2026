#!/usr/bin/env Rscript
#==================================================================================
# Step 10: Genomic Distribution Analysis (Peak Annotation)
# Annotates peaks to genomic features using ChIPseeker
# Tests for distribution changes across treatment conditions
#
# Input:  Peak files from Step 04 (narrowPeak format)
# Output: Distribution tables, statistical tests, visualization PDFs
#
# Dependencies: ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db
#==================================================================================

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(tibble)
})

# Setup
projPath <- Sys.getenv("PROJ_DIR", "/path/to/your/project")
peakDir  <- file.path(projPath, "peakCalling/MACS2/human/Untreated_control")
outDir   <- file.path(projPath, "Analysis/Genomic_Distribution_Detailed")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

conditions <- c("D7_Cis", "D14_Cis", "D14_C70")
condition_labels <- c("Day 7", "Day 14", "Day 14 + C70")
marks <- c("H3K4me3", "KDM5A", "KDM5B")
replicates <- c("rep1", "rep2")

# Color schemes
annotation_colors <- c(
  "Promoter (<=1kb)" = "#E41A1C", "Promoter (1-2kb)" = "#FB8072",
  "Promoter (2-3kb)" = "#FDB462", "5' UTR" = "#377EB8",
  "3' UTR" = "#4DAF4A", "1st Exon" = "#984EA3",
  "Other Exon" = "#FF7F00", "1st Intron" = "#FFFF33",
  "Other Intron" = "#A65628", "Downstream (<=3kb)" = "#F781BF",
  "Distal Intergenic" = "#999999"
)

mark_colors <- c("H3K4me3" = "#E41A1C", "KDM5A" = "#4DAF4A", "KDM5B" = "#984EA3")

simple_categories <- c(
  "Promoter" = "#2171B5", "Gene Body" = "#41AB5D",
  "Downstream" = "#F781BF", "Intergenic" = "#999999"
)

# Load peak files
cat("=== LOADING PEAKS ===\n")
peak_files <- list()
peak_counts <- data.frame()

for(condition in conditions) {
  for(mark in marks) {
    all_peaks <- GRanges()
    for(rep in replicates) {
      peak_file <- file.path(peakDir, condition, mark, rep,
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
      merged <- reduce(all_peaks)
      sample_id <- paste(condition, mark, sep="_")
      peak_files[[sample_id]] <- merged
      peak_counts <- rbind(peak_counts, data.frame(
        condition = condition, mark = mark, peak_count = length(merged)))
    }
  }
}

write.csv(peak_counts, file.path(outDir, "Peak_Counts.csv"), row.names = FALSE)

# Annotate peaks
cat("=== ANNOTATING PEAKS ===\n")
anno_list <- list()
for(sample_id in names(peak_files)) {
  anno <- annotatePeak(peak_files[[sample_id]], tssRegion = c(-3000, 1000),
                       TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
  anno_list[[sample_id]] <- anno
  write.csv(as.data.frame(anno),
            file.path(outDir, paste0(sample_id, "_annotation.csv")), row.names = FALSE)
}

# Calculate distributions
cat("=== CALCULATING DISTRIBUTIONS ===\n")
all_dist_detailed <- data.frame()
all_dist_simple <- data.frame()

for(sample_id in names(anno_list)) {
  parts <- strsplit(sample_id, "_")[[1]]
  condition <- paste(parts[1:(length(parts)-1)], collapse="_")
  mark <- parts[length(parts)]
  anno_df <- as.data.frame(anno_list[[sample_id]])

  # Detailed (11 categories)
  detailed <- anno_df %>%
    mutate(feature = case_when(
      grepl("Promoter \\(<=1kb\\)", annotation) ~ "Promoter (<=1kb)",
      grepl("Promoter \\(1-2kb\\)", annotation) ~ "Promoter (1-2kb)",
      grepl("Promoter \\(2-3kb\\)", annotation) ~ "Promoter (2-3kb)",
      grepl("5' UTR", annotation) ~ "5' UTR",
      grepl("3' UTR", annotation) ~ "3' UTR",
      grepl("1st Exon", annotation) ~ "1st Exon",
      grepl("Other Exon", annotation) ~ "Other Exon",
      grepl("1st Intron", annotation) ~ "1st Intron",
      grepl("Other Intron", annotation) ~ "Other Intron",
      grepl("Downstream", annotation) ~ "Downstream (<=3kb)",
      grepl("Distal Intergenic", annotation) ~ "Distal Intergenic",
      TRUE ~ "Other"
    )) %>%
    group_by(feature) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(100 * count / sum(count), 2),
           condition = condition, mark = mark)
  all_dist_detailed <- rbind(all_dist_detailed, detailed)

  # Simple (4 categories)
  simple <- anno_df %>%
    mutate(feature = case_when(
      grepl("Promoter", annotation) ~ "Promoter",
      grepl("Exon|Intron|UTR", annotation) ~ "Gene Body",
      grepl("Downstream", annotation) ~ "Downstream",
      grepl("Distal", annotation) ~ "Intergenic",
      TRUE ~ "Other"
    )) %>%
    group_by(feature) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(100 * count / sum(count), 2),
           condition = condition, mark = mark)
  all_dist_simple <- rbind(all_dist_simple, simple)
}

write.csv(all_dist_detailed, file.path(outDir, "Distribution_Detailed.csv"), row.names = FALSE)
write.csv(all_dist_simple, file.path(outDir, "Distribution_Simple.csv"), row.names = FALSE)

# Statistical testing
cat("=== STATISTICAL TESTING ===\n")
stat_results <- list()
for(mark in marks) {
  mark_data <- all_dist_simple %>%
    filter(mark == !!mark) %>%
    select(condition, feature, count) %>%
    pivot_wider(names_from = condition, values_from = count, values_fill = 0)
  if(nrow(mark_data) > 0 && ncol(mark_data) > 1) {
    count_matrix <- as.matrix(mark_data[, -1])
    rownames(count_matrix) <- mark_data$feature
    chisq_result <- chisq.test(count_matrix)
    stat_results[[mark]] <- list(test = chisq_result, matrix = count_matrix)
    cat(sprintf("  %s: chi2=%.2f, p=%.4e %s\n", mark,
                chisq_result$statistic, chisq_result$p.value,
                ifelse(chisq_result$p.value < 0.05, "***", "")))
  }
}

# Visualizations
cat("=== CREATING VISUALIZATIONS ===\n")

# Per-marker detailed stacked bars
for(mark in marks) {
  mark_data <- all_dist_detailed %>% filter(mark == !!mark)
  if(nrow(mark_data) > 0) {
    p <- ggplot(mark_data, aes(x = factor(condition, levels = conditions),
                               y = percentage, fill = feature)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = 14) +
      labs(title = paste(mark, "Genomic Distribution"), x = "Treatment", y = "Peaks (%)") +
      scale_fill_manual(values = annotation_colors, drop = FALSE) +
      scale_x_discrete(labels = condition_labels)
    ggsave(file.path(outDir, paste0(mark, "_Distribution_Detailed.pdf")), p, width = 12, height = 8)
  }
}

# All markers comparison
p_all <- ggplot(all_dist_simple, aes(x = factor(condition, levels = conditions),
                                     y = percentage, fill = feature)) +
  geom_bar(stat = "identity") + facet_wrap(~ mark, ncol = 3) +
  theme_minimal(base_size = 12) +
  labs(title = "All Markers Comparison", x = "Treatment", y = "Peaks (%)") +
  scale_fill_manual(values = simple_categories) +
  scale_x_discrete(labels = c("D7", "D14", "C70"))
ggsave(file.path(outDir, "All_Markers_Comparison.pdf"), p_all, width = 14, height = 6)

# Heatmap
heatmap_data <- all_dist_detailed %>%
  mutate(sample = paste(mark, condition, sep = "_")) %>%
  select(sample, feature, percentage) %>%
  pivot_wider(names_from = sample, values_from = percentage, values_fill = 0) %>%
  column_to_rownames("feature")

pdf(file.path(outDir, "Distribution_Heatmap.pdf"), width = 12, height = 8)
pheatmap(as.matrix(heatmap_data), cluster_rows = FALSE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.1f",
         color = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(100),
         main = "Genomic Distribution (%) - All Markers")
dev.off()

cat("\nAnalysis complete! Results in:", outDir, "\n")
