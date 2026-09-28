# 2. snRNA-seq normalisation, dimensionality reduction and clustering
# SCTransform, PCA, UMAP, clustering, cluster markers and gene set enrichment.

library(Seurat)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(openxlsx)

setwd("/path/to/Multiome_project/")
options(future.globals.maxSize = 8000 * 1024^2)

nb.filt.seurat <- readRDS("datafiles/nb_filt_seurat_qc.RDS")

## Normalisation ----
#SCTransform on the top 2000 variable genes, regressing out cell cycle and mitochondrial proportion
nb.filt.seurat$Seurat.Cycle.Score <- nb.filt.seurat$Seurat.S - nb.filt.seurat$Seurat.G2M
set.seed(1204)
nb.filt.seurat <- SCTransform(nb.filt.seurat, method = "glmGamPoi",
                              vst.flavor = "v2",
                              vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"),
                              do.scale = TRUE,
                              do.center = TRUE,
                              variable.features.n = 2000,
                              conserve.memory = TRUE)

## PCA, UMAP and clustering ----
nb.filt.seurat <- RunPCA(nb.filt.seurat, npcs = 100, features = VariableFeatures(nb.filt.seurat), verbose = FALSE)
nb.filt.seurat <- RunUMAP(nb.filt.seurat, dims = 1:30)
nb.filt.seurat <- FindNeighbors(nb.filt.seurat, dims = 1:30)
nb.filt.seurat <- FindClusters(nb.filt.seurat, resolution = 0.2)

DimPlot(nb.filt.seurat, group.by = "seurat_clusters", label = TRUE)
ggsave("plots/rna_umap_clusters.pdf", width = 7, height = 7)
DimPlot(nb.filt.seurat, group.by = "Sample")
ggsave("plots/rna_umap_samples.pdf", width = 7, height = 7)

## Cluster markers ----
Idents(nb.filt.seurat) <- nb.filt.seurat$seurat_clusters
nb.clust.wilcox.markers <- FindAllMarkers(nb.filt.seurat,
                                          assay = "SCT",
                                          slot = "data",
                                          only.pos = TRUE,
                                          min.pct = 0.05,
                                          logfc.threshold = 0.2,
                                          test.use = "wilcox")

nb.clust.wilcox.list <- lapply(unique(nb.clust.wilcox.markers$cluster), function(x){
  tmp.df <- nb.clust.wilcox.markers[nb.clust.wilcox.markers$cluster == x & nb.clust.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE), ]
})
names(nb.clust.wilcox.list) <- unique(nb.clust.wilcox.markers$cluster)

## Gene set enrichment of cluster markers ----
#GO biological process
nb.clust.GO.BP.list <- lapply(nb.clust.wilcox.list, function(x){
  as.data.frame(enrichGO(x$gene, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", readable = TRUE, pAdjustMethod = "BH"))
})

#MSigDB Hallmarks
H.df <- as.data.frame(msigdbr(species = "Homo sapiens", category = "H")[, c("gs_name", "gene_symbol")])
nb.clust.H.list <- lapply(nb.clust.wilcox.list, function(x){
  as.data.frame(enricher(x$gene, pAdjustMethod = "BH", TERM2GENE = H.df))
})

## Save data files ----
write.xlsx(nb.clust.wilcox.list, "datafiles/nb_seurat_markers_wilcox.xlsx")
write.xlsx(nb.clust.GO.BP.list, "datafiles/nb_seurat_markers_BP_ontology.xlsx")
write.xlsx(nb.clust.H.list, "datafiles/nb_seurat_markers_Hallmarks.xlsx")
saveRDS(nb.filt.seurat, "datafiles/nb_seurat_clustered.RDS")
