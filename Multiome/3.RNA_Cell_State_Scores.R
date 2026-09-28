# 3. snRNA-seq cell state scores
# ADRN and MES signatures (van Groningen et al. 2017) and per-cell AMT state.

library(Seurat)
library(readxl)
library(dplyr)
library(ggplot2)

setwd("/path/to/Multiome_project/")

nb.filt.seurat <- readRDS("datafiles/nb_seurat_clustered.RDS")

## ADRN and MES signatures ----
vanGroningen.df <- as.data.frame(read_xlsx(path = "datafiles/vanGroningen_2017.xlsx", col_names = FALSE))
colnames(vanGroningen.df) <- c("Gene", "Signature")
vanGroningen.sig <- split(vanGroningen.df$Gene, vanGroningen.df$Signature)

nb.filt.seurat <- AddModuleScore(nb.filt.seurat, features = vanGroningen.sig, assay = "SCT", seed = 12345,
                                 name = paste0(names(vanGroningen.sig), ".Sig"))
colnames(nb.filt.seurat@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(nb.filt.seurat@meta.data))

## AMT state ----
#For each cell, t-test of MES versus ADRN gene expression
nb.counts <- GetAssayData(nb.filt.seurat, assay = "SCT", layer = "data")
nb.counts.mes <- nb.counts[rownames(nb.counts) %in% vanGroningen.sig$MES, ]
nb.counts.adrn <- nb.counts[rownames(nb.counts) %in% vanGroningen.sig$ADRN, ]

t.stats <- matrix(NA, nrow = ncol(nb.counts), ncol = 2, dimnames = list(colnames(nb.counts), c("t", "p.value")))
for(cell in 1:ncol(nb.counts)) {
  result <- t.test(nb.counts.mes[, cell], nb.counts.adrn[, cell])
  t.stats[cell, ] <- c(result$statistic, result$p.value)
}

amt.states <- as.data.frame(t.stats) %>%
  mutate(AMT.state = case_when(t > 0 & p.value < 0.05 ~ "MES",
                               t < 0 & p.value < 0.05 ~ "ADRN",
                               TRUE ~ "intermediate"))

nb.filt.seurat$AMT.score <- amt.states$t
nb.filt.seurat$AMT.state <- amt.states$AMT.state
table(nb.filt.seurat$Sample, nb.filt.seurat$AMT.state)

amt.cols <- c("ADRN" = "#990099", "intermediate" = "lightgrey", "MES" = "#F37735")
DimPlot(nb.filt.seurat, group.by = "AMT.state", cols = amt.cols)
ggsave("plots/AMT_state_umap.pdf", width = 7, height = 7)

## Save data files ----
saveRDS(nb.filt.seurat, "datafiles/nb_seurat_AMT.rds")
