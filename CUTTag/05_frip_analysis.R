#!/usr/bin/env Rscript
#==================================================================================
# Step 05: Peak Quality Analysis (FRiP, peak numbers, peak widths)
# Calculates Fraction of Reads in Peaks and quality metrics
#
# Input:  Peak files from Step 04, deduplicated BAMs, QC summary
# Output: Peak quality plots (PDF), summary tables
#==================================================================================

library(tidyverse)
library(GenomicRanges)
library(Rsamtools)
library(ggpubr)

# Set paths - edit or set PROJ_DIR environment variable
projPath <- Sys.getenv("PROJ_DIR", "/path/to/your/project")
peakPath <- file.path(projPath, "peakCalling/MACS2")
bamPath  <- file.path(projPath, "bam_files_dedup")

# Read peak summary from Step 04
peak_summary <- read.table(file.path(peakPath, "peak_summary.txt"),
                           header = TRUE, sep = "\t")

# Filter for default results
peak_summary_filtered <- peak_summary %>%
  filter(Q_value == 0.1, Mode == "default") %>%
  distinct(Replicate, Comparison, Target, .keep_all = TRUE)

target_colors <- c("H3K4me3"="#E41A1C", "KDM5A"="#4DAF4A", "KDM5B"="#984EA3")

# =====================================================
# 1. PEAK NUMBER ANALYSIS
# =====================================================
cat("=== PEAK NUMBER ANALYSIS ===\n")

peak_number_summary <- peak_summary_filtered %>%
  group_by(Target, Comparison) %>%
  summarise(
    Mean_Peaks = mean(Peak_Count),
    SD_Peaks = sd(Peak_Count),
    CV = (SD_Peaks / Mean_Peaks) * 100,
    .groups = "drop"
  ) %>%
  mutate(Treatment = sub("_vs_.*", "", Comparison))

p_peak_number <- peak_summary_filtered %>%
  mutate(Treatment = factor(sub("_vs_.*", "", Comparison),
                            levels = c("D7_Cis", "D14_Cis", "D14_C70"))) %>%
  ggplot(aes(x = Target, y = Peak_Count, fill = Target)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_point(aes(shape = as.factor(Replicate)),
             position = position_dodge(width = 0.3), size = 3) +
  facet_wrap(~Treatment, scales = "free_y") +
  scale_fill_manual(values = target_colors) +
  scale_shape_manual(values = c(16, 17), name = "Replicate") +
  theme_bw(base_size = 10) +
  labs(title = "Peak Numbers by Treatment and Target",
       y = "Number of Peaks", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave("peak_numbers_analysis.pdf", p_peak_number, width = 10, height = 6)

# =====================================================
# 2. PEAK WIDTH ANALYSIS
# =====================================================
cat("=== PEAK WIDTH ANALYSIS ===\n")

get_peak_widths <- function(condition, target, rep) {
  peak_file <- file.path(peakPath, "human/Untreated_control", condition,
                         target, paste0("rep", rep),
                         "macs2_peak_q0.1_peaks.narrowPeak.bed")
  if (file.exists(peak_file)) {
    peaks <- read.table(peak_file, header = FALSE, sep = "\t")
    data.frame(Condition = condition, Target = target, Replicate = rep,
               Width = abs(peaks$V3 - peaks$V2))
  } else NULL
}

all_widths <- list()
for (i in 1:nrow(peak_summary_filtered)) {
  row <- peak_summary_filtered[i, ]
  condition <- strsplit(as.character(row$Comparison), "_vs_")[[1]][1]
  rep_num <- as.numeric(gsub("Rep", "", row$Replicate))
  w <- get_peak_widths(condition, row$Target, rep_num)
  if (!is.null(w)) all_widths[[length(all_widths) + 1]] <- w
}

if (length(all_widths) > 0) {
  peak_width_data <- bind_rows(all_widths)

  p_peak_width <- peak_width_data %>%
    ggplot(aes(x = Target, y = Width, fill = Target)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    facet_wrap(~Condition) +
    scale_y_log10(breaks = c(100, 500, 1000, 5000, 10000),
                  labels = c("100", "500", "1K", "5K", "10K")) +
    scale_fill_manual(values = target_colors) +
    theme_bw(base_size = 10) +
    labs(title = "Peak Width Distribution",
         y = "Peak Width (bp, log scale)", x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  ggsave("peak_width_distribution.pdf", p_peak_width, width = 10, height = 6)
}

# =====================================================
# 3. FRiP CALCULATION
# =====================================================
cat("=== FRiP CALCULATION ===\n")

qc_data <- read.table(file.path(projPath, "QC_final/qc_summary.txt"),
                       header = TRUE, sep = "\t")

calculate_frip <- function(peak_file, bam_file, total_mapped_reads) {
  if (!file.exists(peak_file) || !file.exists(bam_file)) return(NA)

  peaks <- read.table(peak_file, header = FALSE, sep = "\t")
  peak_gr <- GRanges(seqnames = peaks$V1,
                     ranges = IRanges(start = peaks$V2, end = peaks$V3))

  param <- ScanBamParam(which = peak_gr,
                        flag = scanBamFlag(isUnmappedQuery = FALSE))
  reads_in_peaks <- sum(countBam(BamFile(bam_file), param = param)$records)

  (reads_in_peaks / total_mapped_reads) * 100
}

frip_results <- peak_summary_filtered %>%
  rowwise() %>%
  mutate(
    Sample_ID = strsplit(Sample_IDs, "_vs_")[[1]][1],
    Condition = strsplit(Comparison, "_vs_")[[1]][1],
    Rep_Num = as.numeric(gsub("Rep", "", Replicate)),
    Peak_File = file.path(peakPath, "human/Untreated_control",
                          Condition, Target, paste0("rep", Rep_Num),
                          "macs2_peak_q0.1_peaks.narrowPeak.bed"),
    BAM_File = file.path(bamPath, paste0(Sample_ID, "_dedup.bam"))
  ) %>%
  left_join(qc_data %>% select(Sample, Mapped_Reads),
            by = c("Sample_ID" = "Sample")) %>%
  rowwise() %>%
  mutate(FRiP = calculate_frip(Peak_File, BAM_File, Mapped_Reads))

p_frip <- frip_results %>%
  ggplot(aes(x = Target, y = FRiP, fill = Target)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +
  geom_point(aes(shape = as.factor(Replicate)),
             position = position_dodge(width = 0.3), size = 3) +
  facet_wrap(~Condition) +
  scale_fill_manual(values = target_colors) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "orange", alpha = 0.5) +
  theme_bw(base_size = 10) +
  labs(title = "Fraction of Reads in Peaks (FRiP)", y = "FRiP (%)", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave("frip_analysis.pdf", p_frip, width = 10, height = 6)

# Combined figure
combined_plot <- ggarrange(p_peak_number, p_frip, ncol = 1, nrow = 2,
                           common.legend = TRUE, legend = "bottom")
ggsave("peak_analysis_combined.pdf", combined_plot, width = 12, height = 10)

# Save summaries
write.table(peak_number_summary, "peak_number_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Analysis Complete ===\n")
