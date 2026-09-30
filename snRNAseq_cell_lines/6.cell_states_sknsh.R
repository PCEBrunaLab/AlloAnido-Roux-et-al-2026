# 6. ADRN / MES cell states, SK-N-SH
# van Groningen et al. 2017 signatures. Per cell, SCT expression of MES genes
# is compared with ADRN genes (t-test): t > 0 and p < 0.05 = MES,
# t < 0 and p < 0.05 = ADRN, otherwise intermediate.
# Usage: Rscript 6.cell_states_sknsh.R <seurat.rds> <vanGroningen_2017.xlsx> <output.rds>

library(Seurat)
library(readxl)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
seu <- readRDS(args[1]); sig_file <- args[2]; out_file <- args[3]

## Signatures ----
vanGroningen.df <- as.data.frame(read_xlsx(sig_file, col_names = FALSE))
colnames(vanGroningen.df) <- c("Gene", "Signature")
vanGroningen.sig <- split(vanGroningen.df$Gene, vanGroningen.df$Signature)[c("MES", "ADRN")]

seu <- AddModuleScore(seu, features = vanGroningen.sig, assay = "SCT", seed = 12345,
                      name = c("MES.Sig", "ADRN.Sig"))
colnames(seu@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(seu@meta.data))
seu$AMT.Sig <- seu$MES.Sig - seu$ADRN.Sig

## AMT state: per-nucleus t-test ----
expr <- GetAssayData(seu, assay = "SCT", layer = "data")
mes  <- expr[rownames(expr) %in% vanGroningen.sig$MES, , drop = FALSE]
adrn <- expr[rownames(expr) %in% vanGroningen.sig$ADRN, , drop = FALSE]

t.stats <- matrix(NA, nrow = ncol(expr), ncol = 2, dimnames = list(colnames(expr), c("t", "p.value")))
for (cell in 1:ncol(expr)) {
  result <- t.test(mes[, cell], adrn[, cell])
  t.stats[cell, ] <- c(result$statistic, result$p.value)
}

amt <- as.data.frame(t.stats) %>%
  mutate(AMT.state = case_when(t > 0 & p.value < 0.05 ~ "MES",
                               t < 0 & p.value < 0.05 ~ "ADRN",
                               TRUE ~ "intermediate"))
seu$AMT.score <- amt$t
seu$AMT.state <- amt$AMT.state
table(seu$Condition, seu$AMT.state)

DimPlot(seu, group.by = "AMT.state", raster = TRUE,
        cols = c(ADRN = "#990099", intermediate = "lightgrey", MES = "#F37735"))
ggsave("umap_AMT_state.pdf", width = 7, height = 6)

saveRDS(seu, out_file)
