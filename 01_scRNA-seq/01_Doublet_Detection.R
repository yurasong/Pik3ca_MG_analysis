#####################################################################################
# Doublet detection: Used for QC step
#
# Description:
#   Detects and removes doublets from scRNA-seq Seurat object using two methods:
#     1. scds (cxds, bcds, and hybrid scoring)
#     2. DoubletFinder
#   Retains cells classified as singlets by both tools, then re-clusters and projects UMAP.
#   This step is done on individual data level, not after integration.
#
# Inputs:
#   - seurat_object.rds       : Seurat object after filtering empty droplets and low-expressed cells
#
# Outputs:
#   - Violin plots of doublet scores (hybrid_score and pANN)
#   - UMAP plot of filtered cells (PDF via DimPlot)
#   - Annotated Seurat object without doublets: Doublet_filtered_object.rds
#
# Dependencies:
#   Seurat, Signac, scds, SingleCellExperiment, magrittr,
#   data.table, dplyr, ggplot2, DoubletFinder, forcats
#
# Reference:
#   - scds: https://github.com/compbiomed/scds
#   - DoubletFinder tutorial: https://github.com/chris-mcginnis-ucsf/DoubletFinder
#####################################################################################

# Library

library(Seurat)
library(scds)
library(DoubletFinder)
library(SingleCellExperiment)
library(magrittr)
library(data.table)
library(dplyr)
library(ggplot2)

# Data

seuset <- readRDS("seurat_object.rds")

seuset@active.ident =factor(as.character(seuset@meta.data$seurat_clusters))
names(seuset@active.ident) = rownames(seuset@meta.data)

DefaultAssay(seuset) <- "RNA"

# Run scds
## Convert to sce object

integrated_sce <- as.SingleCellExperiment(seuset)

## Apply scds

integrated_sce <- cxds(integrated_sce)
integrated_sce <- bcds(integrated_sce)
integrated_sce = cxds_bcds_hybrid(integrated_sce)

col_cd <- colData(integrated_sce)
head(cbind(col_cd$cxds_score,col_cd$bcds_score, col_cd$hybrid_score))

integrated_score <- cbind(col_cd$cxds_score,col_cd$bcds_score, col_cd$hybrid_score)
colnames(integrated_score) <- c("coexp_score", "binary_score", "hybrid_score")
head(integrated_score)

integrated_meta <- seuset@meta.data

mer <- merge(integrated_meta, integrated_score, by="row.names")
rownames(mer) <- mer$Row.names
mer <- mer[, -1]

## Visualisation

seuset@meta.data <- mer

fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'blue'),  limits = c(0, 2))
FeaturePlot(seuset, features="hybrid_score", label=F, pt.size = 0.5) +
  ggtitle("scds doublet detection") + fix.sc

# Run DoubletFinder

## Apply DoubletFinder

sweep.res.list <- paramSweep(seuset, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

bcmvn <- find.pK(sweep.stats)

## Homotypic estimation

annotations <- seuset@meta.data$seurat_clusters

nExp_poi <- round(0.06*nrow(seuset@meta.data)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder

seuset <- doubletFinder(seuset, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

## Visualisation

fix.sc <- scale_color_gradientn(colours = c('lightgrey', 'blue'))
FeaturePlot(seuset, features="pANN_0.25_0.09_*", pt.size=0.5) + fix.sc # Change * to your result 