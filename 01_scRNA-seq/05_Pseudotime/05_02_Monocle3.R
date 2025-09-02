##############################################################################################
# 05_02_Monocle3.R
# This code is used for inferring pseudotime trajectory on scRNA-seq data with monocle3.
##############################################################################################

# Preparation
## Package install
### Dependencies

# comment: Those dependencies are not installed automatically.

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

## monocle installation
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

## Library
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(monocle3)

## Data
seurat_obj <- readRDS("../00_Pseudotime_input.rds")

# Apply monocle
## Convert monocle CDS
### Extract expression matrix
data <- as(as.matrix(GetAssayData(seurat_obj, slot = "counts")), "sparseMatrix")

### Extract cell metadata and feature metadata
cell_meta <- as.data.frame(seurat_obj@meta.data[, c(1:4, 6, 8)])
gene_meta <- as.data.frame(rownames(data), row.names=rownames(data))
names(gene_meta) <- "gene_short_name"

## Create Monocle CellDataSet
cds <- new_cell_data_set(data,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 40, method="PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')

#plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "cell_type")


# Pull-out UMAP coordinate
# comment: Use same UMAP coordinate of input data
seurat_umap <- Embeddings(seurat_obj, "umap")
all(rownames(seurat_umap) %in% colnames(cds))  # should be TRUE

# assign them into the 'reducedDims' slot of the cds
reducedDims(cds)[["UMAP"]] <- seurat_umap[ colnames(cds), ]

plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by      = "cell_type",
           label_cell_groups   = FALSE)

cds <- cluster_cells(cds)
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by      = "cell_type",
           label_cell_groups   = FALSE, 
           cell_size=1,
           cell_stroke = 0) + NoLegend() +
  scale_color_manual(values=c("#32CD32", "#800000", "#F08080", "#BDB76B",  "#1E90FF"))


## Fit trajectory
cds <- learn_graph(cds)
plot_cells(cds,
           reduction_method = "UMAP",
           color_cells_by      = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           cell_size=1,
           cell_stroke = 0) +
  scale_color_manual(values=c("#32CD32", "#800000", "#F08080", "#BDB76B",  "#1E90FF"))

# If selection of start and end point is needed, 
cds <- order_cells(cds,  reduction_method = "UMAP") # Will show shiny window.

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

