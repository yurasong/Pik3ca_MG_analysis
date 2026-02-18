#####################################################################################
# SCTransform-based integration of scRNA-seq datasets (CTL, ER_Pik, Kit_Pik)
#
# Description:
#   Integrates three scRNA-seq Seurat objects (CTL, ER_Pik, and Kit_Pik) using
#   SCTransform normalization with regression of mitochondrial content and
#   cell-cycle scores. The pipeline performs cell-cycle scoring per dataset,
#   applies SCT with covariate regression, identifies integration features,
#   computes anchors using RPCA, integrates the datasets, and generates PCA/UMAP
#   embeddings with graph-based clustering. Finally, the RNA layer is re-joined
#   to enable downstream visualization and marker analysis in the RNA assay.
#
# Workflow overview:
#     1. Load individual annotated/filtered Seurat objects (CTL, ER_Pik, Kit_Pik)
#     2. Add sample labels
#     3. Compute S and G2M scores (CellCycleScoring) per dataset
#     4. SCTransform with regression (percent.mt, S.Score, G2M.Score)
#     5. Select integration features and prepare SCT integration
#     6. Find integration anchors (RPCA) and integrate datasets (SCT)
#     7. Run PCA/UMAP, construct neighbors, and perform clustering
#     8. Switch back to RNA assay and JoinLayers for downstream RNA-based plots
#     9. Export integrated object
#
# Inputs:
#   - Control_annotated.rds         : Control (CTL) Seurat object
#   - ERPik_annotated.rds           : ER_Pik Seurat object
#   - Kit-Pik_scRNA_filtered.rds    : Kit_Pik Seurat object (filtered)
#
# Outputs:
#   - Seurat_integrated_res0p5.rds
#       Integrated Seurat object (SCT integration) with PCA/UMAP and clusters
#
#
# Dependencies:
#   Seurat, ggplot2, patchwork, tidyverse, data.table, magrittr
#
# Notes:
#   - SCT integration is used to mitigate batch effects while controlling
#     for cell-cycle-driven variance.
#   - JoinLayers is run after integration to ensure RNA-layer availability
#     for downstream differential expression and FeaturePlot/DotPlot.
#####################################################################################


# Library
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

# Individual data
ctl <- readRDS("Control_annotated.rds")
er <- readRDS("ERPik_annotated.rds")
kit <- readRDS("Kit-Pik_scRNA_filtered.rds")

ctl$sample <- 'CTL'
er$sample  <- 'ER_Pik'
kit$sample <- 'Kit_Pik'

# Cell cycle information
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)

ctl <- CellCycleScoring(ctl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
er <- CellCycleScoring(er, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
kit <- CellCycleScoring(kit, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Integration
mar.list <- list(ctl, kit, er)

mar.list <- lapply(mar.list, function(x) {
  SCTransform(x,
              vars.to.regress = c('percent.mt','S.Score','G2M.Score'),
              verbose = FALSE
  )
})

features <- SelectIntegrationFeatures(mar.list, nfeatures = 2000)
mar.list <- PrepSCTIntegration(mar.list, anchor.features = features)

mar.anchors <- FindIntegrationAnchors(
  object.list           = mar.list,
  normalization.method  = 'SCT',
  anchor.features       = features,
  reduction             = 'rpca',
  k.filter              = 30
)

mar.combined <- IntegrateData(
  anchorset            = mar.anchors,
  normalization.method = 'SCT',
  dims                 = 1:30
)

# Clustering 
DefaultAssay(mar.combined) <- 'integrated'
mar.combined <- ScaleData(mar.combined, verbose = FALSE)
mar.combined <- RunPCA(mar.combined, npcs = 30, verbose = FALSE)
mar.combined <- RunUMAP(mar.combined, reduction = 'pca', dims = 1:30)
mar.combined <- FindNeighbors(mar.combined, reduction = 'pca', dims = 1:30)

mar.combined_res0p5 <- FindClusters(mar.combined, resolution = 0.5)

# Join RNA layer
DefaultAssay(mar.combined_res0p5) <- "RNA"
mar.combined_res0p5 <- JoinLayers(mar.combined_res0p5)

# Export result
saveRDS(mar.combined_res0p5, "Seurat_integrated_res0p5.rds")