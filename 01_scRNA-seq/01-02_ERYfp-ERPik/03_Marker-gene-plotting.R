#####################################################################################
# UMAP feature plots for stromal marker genes (Extended Data Fig. 4c–h)
#
# Description:
#   Generates UMAP-based expression overlays for selected stromal lineage marker
#   genes using the integrated and annotated stromal Seurat object. The script
#   extracts UMAP coordinates and log-normalized RNA expression values, then
#   produces ggplot-based feature maps with a consistent viridis colour scale
#   and minimalist theme for figure-panel assembly.
#
# Workflow overview:
#     1. Load integrated annotated Seurat object (stromal + epithelial populations)
#     2. Define a viridis continuous colour palette for expression values
#     3. Extract UMAP embeddings and merge with RNA assay expression matrix
#     4. Define helper function to plot one gene on UMAP (custom feature plot)
#     5. Plot marker gene sets for major stromal cell types (Extended Data Fig. 4c–h)
#
# Inputs:
#   - Stroma_Integrated_Annot.rds
#       Integrated Seurat object with UMAP embedding and `cell_type` annotations
#
# Outputs:
#   - A list of ggplot objects (one per gene) showing expression on UMAP
#     (typically saved/arranged downstream for Extended Data Fig. 4c–h)
#
# Marker genes (figure panels):
#   - Fibroblast : Col1a1, Pdgfra
#   - B cell     : Cd79a, Fcmr
#   - T cell     : Cd3g, Icos
#   - Macrophage : Cd86, Csf1r
#   - Endothelial: Pecam1, Cldn5
#   - Pericyte   : Cspg4, Myh11
#
# Dependencies:
#   Seurat, tidyverse, dplyr, data.table, ggplot2, viridis,
#   RColorBrewer, plotly, grid
#
# Notes:
#   - Expression is taken from `seuset[["RNA"]]$data` (log-normalized values).
#   - Legends and axis elements are removed to facilitate figure assembly.
#   - This custom plotting approach mirrors FeaturePlot but ensures consistent
#     aesthetics across panels.
#####################################################################################


# Library

library(Seurat)

library(RColorBrewer)
library(dplyr)
library(plotly)
library(grid)
library(data.table)
library(tidyverse)

library(viridis)

# Data

seuset <- readRDS("Stroma_Integrated_Annot.rds")

# Plotting
## Colour palette

myPalette <- viridis(n = 10, option = "C", direction = -1)
sc <- scale_colour_gradientn(colours = rev(myPalette))

## Set-up

umap_coord <- as.data.frame(seuset@reductions$umap@cell.embeddings)
seuset_gene <- merge.data.frame(umap_coord,t(seuset[["RNA"]]$data),by=0)
row.names(seuset_gene) <- seuset_gene$Row.names

## Code for visualisation

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- seuset_gene[,regulon]
  regulon_plot <- ggplot(data=seuset_gene) + geom_point(aes(x=umap_1,y=umap_2,colour=regulon_d), size=1)
  regulon_plot <- regulon_plot + labs(colour=regulon)
  regulon_plot <- regulon_plot + sc
  regulon_plot <- regulon_plot + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  regulon_plot
  
}


# Marker genes shown on Extended Fig. 4c-h

fibro <- c("Col1a1", "Pdgfra")
bcell <- c("Cd79a", "Fcmr")
tcell <- c("Cd3g", "Icos")
macro <- c("Cd86", "Csf1r")
endo <- c("Pecam1", "Cldn5")
pericyte <- c("Cspg4", "Myh11")

all_gene <- c(bcell, tcell, fibro, macro, endo, pericyte)

lapply(all_gene, umap_color_scaled)
