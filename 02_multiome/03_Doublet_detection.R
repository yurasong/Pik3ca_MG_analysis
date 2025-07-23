#####################################################################################
# File: 03_Doublet_detection.R
#
# Description:
#   Detects and removes doublets from a filtered multiome Seurat object using two methods:
#     1. scds (cxds, bcds, and hybrid scoring)
#     2. DoubletFinder
#   Retains cells classified as singlets by both tools, then re-clusters and projects UMAP.
#
# Inputs:
#   - Annot_multiome_noFibro.rds       : Seurat object after fibroblast filtering
#
# Outputs:
#   - Violin plots of doublet scores (hybrid_score and pANN)
#   - UMAP plot of filtered cells (PDF via DimPlot)
#   - Annotated Seurat object without doublets: Annot_multiome_noFibro_nodoublet.rds
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
library(Signac)
library(scds)
library(SingleCellExperiment)
library(magrittr)
library(data.table)
library(dplyr)
library(ggplot2)
library(DoubletFinder)

# Data
seuset <- readRDS("Annot_multiome_noFibro.rds")
seuset@active.ident =factor(as.character(seuset@meta.data$seurat_clusters))
names(seuset@active.ident) = rownames(seuset@meta.data)

DefaultAssay(seuset) <- "RNA"

# Running scds
## Need to convert to SingleCellExperiment data from seurat object.
integrated_sce <- as.SingleCellExperiment(seuset)

## Applying methodology
integrated_sce <- cxds(integrated_sce)
integrated_sce <- bcds(integrated_sce)

integrated_sce = cxds_bcds_hybrid(integrated_sce)

col_cd <- colData(integrated_sce)

## Put back to seurat object
integrated_score <- cbind(col_cd$cxds_score,col_cd$bcds_score, col_cd$hybrid_score)
colnames(integrated_score) <- c("coexp_score", "binary_score", "hybrid_score")
head(integrated_score)

integrated_meta <- seuset@meta.data

mer <- merge(integrated_meta, integrated_score, by="row.names")
rownames(mer) <- mer$Row.names
mer <- mer[, -1]

seuset@meta.data <- mer

## Plotting

VlnPlot(seuset, features="hybrid_score", pt.size=0) + 
  ggtitle("scds doublet detection on Integrated") + 
  NoLegend() +
  ylim(0, 2) +
   geom_hline(yintercept=0.75, linetype="dashed", color = "black") # threshold of doublet is 0.75.


# Running DoubletFinder

sweep.res.list_cd1 <- paramSweep_v3(seuset, PCs = 1:10, sct = FALSE)
sweep.stats_cd1 <- summarizeSweep(sweep.res.list_cd1, GT = FALSE)
bcmvn_cd1 <- find.pK(sweep.stats_cd1)

annotations <- seuset@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_cd1@meta.data$ClusteringResults
nExp_poi <- round(0.06*nrow(seuset@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Add score and classification
seuset <- doubletFinder_v3(seuset, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

head(seuset@meta.data) # to validate

VlnPlot(seuset, features="pANN_0.25_0.09_664", pt.size=0) + 
  ggtitle("doubletFinder on CTL") + 
  NoLegend() +
  ylim(0, 0.7) +
  geom_hline(yintercept=0.4, linetype="dashed", color = "black") # score starts with pANN_, but name could be different per run.

# Filter out the scores
seuset_filt <- seuset[, seuset$pANN_0.25_0.09_664 < 0.5 & seuset$hybrid_score < 0.75]

## Re-clustering
### Gene expression

DefaultAssay(seuset_filt) <- "RNA"
subc <- FindVariableFeatures(seuset_filt, selection.method = "vst", nfeatures = 1000)

all.genes <- rownames(subc)
subc <- ScaleData(subc, features = all.genes)

subc <- RunPCA(subc, features = VariableFeatures(object = subc))
subc <- FindNeighbors(subc, dims = 1:40, k.param = 5)

### DNA accessibility data processing

DefaultAssay(subc) <- "peaks"
subc <- FindTopFeatures(subc, min.cutoff = 5)
subc <- RunTFIDF(subc)
subc <- RunSVD(subc) 

### Joint UMAP visualization
subc <- FindMultiModalNeighbors(
  object = subc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:40), 
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

DefaultAssay(subc) <- "RNA"
subc <- FindClusters(subc, resolution = 0.4)

subc <- RunUMAP(
  object = subc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(subc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend() 

# Annotation

subc <- RenameIdents(subc, 
                       `0` = "LC_ER-", 
                       `1` = "LC_ER+", 
                       `2` = "LC_ER+",
                       `3` = "LC_ER-", 
                       `4` = "LC_ER-", 
                       `5` = "BCs")

subc$cell_type <- fct_infreq(subc@active.ident)

subc@active.ident =factor(as.character(subc@meta.data$cell_type))
names(subc@active.ident) = rownames(subc@meta.data)

DimPlot(subc, reduction = "umap", label = TRUE, label.size = 4, repel = T,
        cols=c("LC_ER-" = "#00BFC4", "BCs" = "#C77CFF", "LC_ER+" = "#00A9FF", "Prolif" = "#FF61CC")) + 
  NoLegend() 

saveRDS(subc, "Annot_multiome_noFibro_nodoublet.rds")






