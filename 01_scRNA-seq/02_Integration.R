#####################################################################################
# 03_Integration.R
# This code is used for integration of 10X scRNA-seq datasets.
#####################################################################################

# Library
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

# Individual data
ctl <- readRDS("../rds/Control_annotated.rds")
er <- readRDS("../rds/ERPik_annotated.rds")
kit <- readRDS("../rds/Kit-Pik_scRNA_filtered.rds")

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