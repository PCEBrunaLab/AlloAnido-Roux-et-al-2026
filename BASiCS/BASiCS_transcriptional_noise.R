# ============================================================
# Title:        Transcriptional noise quantification with BASiCS
# Input:        pre-processed scRNA-seq object (Seurat or SingleCellExperiment)
# Output:       chains/          BASiCS MCMC chains
#               results/         differential over-dispersion and HVG tables, plots
# Depends on:   R >= 4.2, BASiCS, SingleCellExperiment, Seurat, dplyr, ggplot2
# Notes:        Example shown for MES cells across the cisplatin course.
#               The same steps apply to comparisons between cell states in
#               untreated cells (change GROUPS in section 1).
# ============================================================

library(Seurat)
library(SingleCellExperiment)
library(BASiCS)
library(dplyr)
library(ggplot2)

set.seed(100)

INPUT_RDS   <- "Roux_SKNSH.rds"     # pre-processed object (see README)
CHAIN_DIR   <- "chains"
RESULTS_DIR <- "results"
N_CELLS     <- 130                  # cells per group: size of the smallest group
dir.create(CHAIN_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)


## 1. Data preparation ------------------------------------------------------

seu <- readRDS(INPUT_RDS)
sce <- as.SingleCellExperiment(seu)

COURSE <- c("Untreated", "Cisplatin(1)_ON", "Cisplatin(2)_ON",
            "Cisplatin_1weeksOFF", "Cisplatin_4weeksOFF")
sce_cis <- sce[, sce$Condition %in% COURSE]

# Replicate (A, B, C) as BASiCS batch.
sce_cis$BatchInfo <- recode(sub(".* ", "", as.character(sce_cis$description)),
                            "A" = "1", "B" = "2", "C" = "3")

# Remove genes with no expression in any cell.
sce_cis <- sce_cis[rowSums(counts(sce_cis)) > 0, ]

# Groups to compare: MES cells at each phase of the course.
GROUPS <- c(ut  = "Untreated",
            cis = "Cisplatin(1)_ON",
            r1  = "Cisplatin_1weeksOFF",
            r4  = "Cisplatin_4weeksOFF")
LABELS <- c(ut = "Untreated", cis = "Cisplatin",
            r1 = "1-week recovery", r4 = "4-week recovery")

sce_groups <- lapply(GROUPS, function(cond) {
    x <- sce_cis[, sce_cis$Condition == cond & sce_cis$AMT.state == "MES"]
    x[, sample(ncol(x), N_CELLS)]          # equal number of cells per group
})


## 2. Parameter estimation ----------------------------------------------------

chains <- lapply(names(sce_groups), function(g) {
    BASiCS_MCMC(sce_groups[[g]], N = 4000, Thin = 10, Burn = 2000,
                WithSpikes = FALSE, Regression = FALSE,
                PrintProgress = TRUE, StoreChains = TRUE,
                StoreDir = CHAIN_DIR, RunName = g)
})
names(chains) <- names(sce_groups)

# Stored chains can be reloaded instead of re-running the MCMC:
# chains <- lapply(setNames(names(GROUPS), names(GROUPS)),
#                  BASiCS_LoadChain, StoreDir = CHAIN_DIR)

# Convergence check for mean (mu) and over-dispersion (delta).
pdf(file.path(RESULTS_DIR, "chain_traces.pdf"))
for (g in names(chains)) {
    plot(chains[[g]], Param = "mu",    Gene = 1, log = "y")
    plot(chains[[g]], Param = "delta", Gene = 1, log = "y")
}
dev.off()


## 3. Differential over-dispersion between groups ------------------------------

COMPARISONS <- list(c("ut", "cis"), c("ut", "r1"), c("ut", "r4"),
                    c("cis", "r1"), c("cis", "r4"))

for (cmp in COMPARISONS) {
    a <- cmp[1]; b <- cmp[2]; tag <- paste0(a, "_vs_", b)
    de <- BASiCS_TestDE(Chain1 = chains[[a]], Chain2 = chains[[b]],
                        GroupLabel1 = LABELS[[a]], GroupLabel2 = LABELS[[b]],
                        EpsilonM = log2(2), EpsilonD = log2(1.5), EpsilonR = 0.41,
                        ProbThresholdM = 0.85, Plot = FALSE)
    write.csv(as.data.frame(de, Parameter = "Disp"),
              file.path(RESULTS_DIR, paste0("DiffDisp_", tag, ".csv")), row.names = FALSE)
    ggsave(file.path(RESULTS_DIR, paste0("DiffDisp_volcano_", tag, ".pdf")),
           BASiCS_PlotDE(de, Plots = "Volcano", Parameters = "Disp"),
           width = 6, height = 5)
}


## 4. Highly variable genes within each group ----------------------------------

for (g in names(chains)) {
    hvg <- BASiCS_DetectHVG(chains[[g]], VarThreshold = 0.6)
    write.csv(as.data.frame(hvg@Table),                      # all genes, HVG TRUE/FALSE
              file.path(RESULTS_DIR, paste0("HVG_all_", g, ".csv")), row.names = FALSE)
    write.csv(as.data.frame(hvg),                            # HVG only
              file.path(RESULTS_DIR, paste0("HVG_", g, ".csv")), row.names = FALSE)
    ggsave(file.path(RESULTS_DIR, paste0("HVG_plot_", g, ".pdf")),
           BASiCS_PlotVG(hvg, "VG"), width = 6, height = 5)
}

sessionInfo()
