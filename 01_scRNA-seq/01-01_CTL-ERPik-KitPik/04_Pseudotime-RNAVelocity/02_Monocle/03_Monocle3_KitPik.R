#####################################################################################
# Monocle3 trajectory inference and pseudotime visualization (Kit Pik Yfp epithelial subset)
#
# Description:
#   Performs trajectory inference using Monocle3 on the integrated epithelial
#   Seurat object. Cells annotated as proliferative ("Prolif") are excluded,
#   the Seurat object is converted to a Monocle3 cell_data_set (CDS) using raw
#   counts, and dimensional reduction (PCA→UMAP) is computed. A Kit Pik Yfp-only subset
#   is then selected to learn a principal graph, identify branch points, and
#   order cells along pseudotime. Outputs are UMAP visualizations colored by
#   cell type and pseudotime for figure panels / sanity checks.
#
# Workflow overview:
#     1. Load integrated epithelial Seurat object
#     2. Subset to non-proliferative cells (exclude "Prolif")
#     3. Convert Seurat counts + metadata into Monocle3 CDS
#     4. Preprocess CDS (PCA) and compute UMAP embedding
#     5. Subset CDS to Kit Pik Yfp cells for trajectory learning
#     6. Cluster cells (Monocle3), learn principal graph, and visualize topology
#     7. Order cells along the learned graph and plot pseudotime on UMAP
#
# Inputs:
#   - Epithelial_Integrated_Annot.rds
#       Integrated epithelial Seurat object with `cell_type` and `sample`
#
# Outputs (Extended Data Fig. 4e):
#   - UMAP plots colored by:
#       • cell_type (pre- and post-graph learning; with/without branch points)
#       • pseudotime (ordered cells along trajectory)
#
# Dependencies:
#   Seurat, monocle3, tidyverse, dplyr, ggplot2
#
# Notes:
#   - CDS is constructed from raw RNA counts (sparseMatrix), not normalized data.
#   - The trajectory is learned on Kit Pik Yfp cells only (sample == "Kit_Pik") to define
#     a baseline manifold; interpretation of pseudotime should reflect this choice.
#   - Ensure `cell_meta` contains the required columns (including `cell_type`
#     and `sample`) and that rownames match cell barcodes.
#####################################################################################


# Library
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(monocle3)

# Data
seuset <- readRDS("Epithelial_Integrated_Annot.rds")
seurat_obj <- seuset[, seuset$cell_type != "Prolif"]

## Convert monocle CDS
data <- as(as.matrix(GetAssayData(seurat_obj, slot = "counts")), "sparseMatrix")
cell_meta <- as.data.frame(seurat_obj@meta.data[, c(1:4, 8, 11)])
gene_meta <- as.data.frame(rownames(data), row.names=rownames(data))
names(gene_meta) <- "gene_short_name"

cds <- new_cell_data_set(data,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)

## Pre-process the data
cds <- preprocess_cds(cds, num_dim = 40, method="PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "cell_type")

# Selection of sample for UMAP calculation
cds_erp <- cds[, colData(cds)$sample == "Kit_Pik"]
plot_cells(cds_erp, label_groups_by_cluster=FALSE, color_cells_by = "cell_type")
cds_erp <- cluster_cells(cds_erp)

## UMAP generation
plot_cells(cds_erp,
           reduction_method = "UMAP",
           color_cells_by      = "cell_type",
           label_cell_groups   = FALSE, 
           cell_size=1,
           cell_stroke = 0) + NoLegend() +
  scale_color_manual(values=c("#FFA500", "#2F4F4F", "#1E90FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## Fit trajectory
cds_erp <- learn_graph(cds_erp)

plot_cells(cds_erp,
           reduction_method = "UMAP",
           color_cells_by      = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

plot_cells(cds_erp,
           reduction_method = "UMAP",
           color_cells_by = "cell_type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           cell_size=1,
           cell_stroke = 0) +
  scale_color_manual(values=c("#FFA500", "#2F4F4F", "#1E90FF")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 

## Set Order (Extended Data Fig. 4e)
cds_erp <- order_cells(cds_erp,  reduction_method = "UMAP")

plot_cells(cds_erp,
           reduction_method = "UMAP",
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
