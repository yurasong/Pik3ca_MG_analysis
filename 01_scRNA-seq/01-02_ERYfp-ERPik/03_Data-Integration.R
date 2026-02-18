#####################################################################################
# Seurat integration and stromal cell-type annotation (ER_Pik and ER_Yfp)
#
# Description:
#   Integrates stromal-containing scRNA-seq datasets from ER Pik Yfp and ER Yfp
#   conditions using Seurat anchor-based integration. After integration, the
#   pipeline performs scaling, PCA, UMAP embedding, neighbor graph construction,
#   clustering, and manual cluster-to-cell-type annotation. UMAPs are generated
#   for (i) global integrated view and (ii) sample-split view (Extended Data Fig. 5a–b).
#
# Workflow overview:
#     1. Load annotated Seurat objects for ER_Yfp (control) and ER_Pik (oncogene)
#     2. Add sample labels and merge into a single object
#     3. Split by sample, normalize, and identify variable features per dataset
#     4. Select integration features and compute integration anchors
#     5. Integrate data, scale, run PCA/UMAP, build neighbors, and cluster
#     6. Rename clusters to biological cell types (manual annotation)
#     7. Plot UMAP (combined + split-by-sample) and export annotated object
#
# Inputs:
#   - ERYFP_Stroma_annotate.rds      : Seurat object (control; stromal included)
#   - ERPik_Stroma_annotate.rds      : Seurat object (ER Pik Yfp; stromal included)
#
# Outputs:
#   - UMAP (integrated; labeled) used for Extended Data Fig. 5a
#   - UMAP (integrated; split by sample) used for Extended Data Fig. 5b
#   - Stroma_Integrated_Annot.rds    : integrated + annotated Seurat object
#
#
# Dependencies:
#   Seurat, ggplot2, patchwork, tidyverse, data.table, magrittr
#
# Notes:
#   - Integration is performed on the "integrated" assay; downstream dimensional
#     reduction and clustering are run in this space.
#   - Cell type labels are assigned manually and stored in `cell_type`
#     (ordered by frequency using fct_infreq).
#####################################################################################


# Library

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

# Data

ctl <- readRDS("ERYFP_Stroma_annotate.rds")
pik <- readRDS("ERPik_Stroma_annotate.rds")

ctl$sample <- "ER_Yfp"
pik$sample <- "ER_Pik"

# Preparation
## Merge dataset

m1 <- merge(pik, ctl, add.cell.ids=c("ER_Pik", "ER_Yfp"))

mar.list <- SplitObject(m1, split.by = "sample")

mar.list <- lapply(X = mar.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = mar.list) # select features that are repeatedly variable across datasets for integration

# Integration

mar.anchors <- FindIntegrationAnchors(object.list = mar.list, anchor.features = features)
mar.combined <- IntegrateData(anchorset = mar.anchors)

DefaultAssay(mar.combined) <- "integrated"
mar.combined <- ScaleData(mar.combined, verbose = FALSE)
mar.combined <- RunPCA(mar.combined, npcs = 40, verbose = FALSE)
mar.combined <- RunUMAP(mar.combined, reduction = "pca", dims = 1:40)
mar.combined <- FindNeighbors(mar.combined, reduction = "pca", dims = 1:40)

mar.combined5 <- FindClusters(mar.combined, resolution = 0.5)

# Annotation

seuset <- RenameIdents(mar.combined5, 
                       `0` = "Fibroblast", 
                       `1` = "BC/Myo", 
                       `2` = "Fibroblast",
                       `3` = "T-cell", 
                       `4` = "B-cell", 
                       `5` = "Fibroblast", 
                       `6` = "LC_ER-",
                       `7` = "Macrophage", 
                       `8` = "HY_BC/ER+",
                       `9` = "Fibroblast",
                       `10` = "Endothelial", 
                       `11` = "LC_ER+",
                       `12` = "T-cell",
                       `13` = "Interferon_producing_Dendritic", 
                       `14` = "Pericyte/SMA",
                       `15` = "Fibroblast",
                       `16` = "Macrophage")

seuset$cell_type <- fct_infreq(seuset@active.ident)

seuset@active.ident =factor(as.character(seuset@meta.data$cell_type))
names(seuset@active.ident) = rownames(seuset@meta.data)

## UMAP used on Extended Data Fig. 5a
DimPlot(seuset, reduction="umap", label=T, repel=T, 
        cols=c("BC/Myo" = "#1E90FF",
               "LC_ER+" = "#FF6347",
               "LC_ER-" = "#2E8B57",
               "Fibroblast" = "#e798f2",
               "B-cell" = "#b7b8fd",
               "T-cell" = "#9a1fae",
               "Macrophage"= "#99e2af",
               "Interferon_producing_Dendritic" = "#8cc5e4",
               "HY_BC/ER+" = "#32CD32",
               "Endothelial" = "#cead2d",
               "Pericyte/SMA" = "#c37397"
        )) + 
  NoLegend()

## UMAP used on Extended Data Fig. 5b
DimPlot(seuset, reduction="umap", label=T, repel=T, split.by="sample",
        cols=c("BC/Myo" = "#1E90FF",
               "LC_ER+" = "#FF6347",
               "LC_ER-" = "#2E8B57",
               "Fibroblast" = "#e798f2",
               "B-cell" = "#b7b8fd",
               "T-cell" = "#9a1fae",
               "Macrophage"= "#99e2af",
               "Interferon_producing_Dendritic" = "#8cc5e4",
               "HY_BC/ER+" = "#32CD32",
               "Endothelial" = "#cead2d",
               "Pericyte/SMA" = "#c37397"
        )) + 
  NoLegend()

saveRDS(seuset, "Stroma_Integrated_Annot.rds")

