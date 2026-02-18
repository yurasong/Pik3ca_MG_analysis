#####################################################################################
# SCTransform-based integration and epithelial state annotation (CTL, ER_Pik, Kit_Pik)
#
# Description:
#   Integrates three epithelial scRNA-seq datasets (CTL, ER_Pik, Kit_Pik) using
#   Seurat SCTransform-based integration while regressing out mitochondrial
#   content and cell-cycle effects. After integration, the pipeline performs
#   dimensional reduction (PCA/UMAP), graph-based clustering, and manual
#   annotation of epithelial states (ER+ LC subsets, ER- LC subsets, HY states,
#   Myo, Prol). UMAPs are generated for the integrated dataset (Fig. 2a) and
#   split by sample/condition (Fig. 2b).
#
# Workflow overview:
#     1. Load individual annotated/filtered Seurat objects (CTL, ER_Pik, Kit_Pik)
#     2. Add sample labels
#     3. Compute cell-cycle scores (S and G2M) per dataset
#     4. SCTransform per dataset with regression (percent.mt, S.Score, G2M.Score)
#     5. Select integration features, prep SCT integration, and find anchors (RPCA)
#     6. Integrate datasets (SCT), run PCA/UMAP, build neighbors, and cluster
#     7. Join RNA layers for downstream RNA-based visualization/DE
#     8. Manually annotate clusters into epithelial states
#     9. Generate UMAPs for figure panels (integrated + split-by-sample)
#
# Inputs:
#   - Control_annotated.rds             : Control (CTL) epithelial Seurat object
#   - ERPik_annotated.rds               : ER_Pik epithelial Seurat object
#   - Kit-Pik_scRNA_filtered.rds        : Kit_Pik epithelial Seurat object (filtered)
#
# Outputs:
#   - UMAP of integrated data (Fig. 2a)
#   - UMAP of integrated data split by sample (Fig. 2b)
#   - Integrated annotated object: Epithelial_Integrated_Annot.rds
#
#
# Dependencies:
#   Seurat, ggplot2, patchwork, tidyverse, data.table, magrittr
#
# Notes:
#   - SCT integration is used to reduce batch effects while controlling for
#     cell-cycle-driven variance.
#   - JoinLayers is run after clustering to ensure RNA-layer availability for
#     downstream DE and FeaturePlot/DotPlot.
#   - Your plotting code suppresses legends/axes for easy multi-panel assembly.
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

# Annotation

seuset <- RenameIdents(mar.combined_res0p5, 
                       `0` = "Late HY", 
                       `1` = "ER+ LCs Sca1", 
                       `2` = "ER- LCs Itgb3low",
                       `3` = "ER- LCs Itgb3high", 
                       `4` = "ER+ LCs Foxa1", 
                       `5` = "ER- LCs Itgb3low", 
                       `6` = "Early HY",
                       `7` = "ER- LCs Itgb3highP",
                       `8` = "ER- LCs Itgb3low",
                       `9` = "Early HY",
                       `10` = "Myo", 
                       `11` = "Prol")

seuset$cell_type <- fct_infreq(seuset@active.ident)

seuset@active.ident =factor(as.character(seuset@meta.data$cell_type))
names(seuset@active.ident) = rownames(seuset@meta.data)

# UMAP plot
## UMAP of integrated data: Fig. 2a

DimPlot(seuset, reduction = "umap", label = F,
        cols = c("Late HY" = "#32CD32", "ER+ LCs Sca1" = "#800000", "ER+ LCs Foxa1" = "#F08080",
                 "ER- LCs Itgb3low" = "#FFA500", "ER- LCs Itgb3high" = "#2F4F4F", 
                 "Early HY" = "#BDB76B", "Myo" = "#1E90FF", "Prol" = "#9932CC")) +
  NoLegend() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## UMAP of integrated data, split.by="sample": Fig. 2b

DimPlot(seuset, reduction = "umap", label = F, split.by="sample",
        cols = c("Late HY" = "#32CD32", "ER+ LCs Sca1" = "#800000", "ER+ LCs Foxa1" = "#F08080",
                 "ER- LCs Itgb3low" = "#FFA500", "ER- LCs Itgb3high" = "#2F4F4F", 
                 "Early HY" = "#BDB76B", "Myo" = "#1E90FF", "Prol" = "#9932CC")) +
  NoLegend() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Export

saveRDS(seuset, "Epithelial_Integrated_Annot.rds")