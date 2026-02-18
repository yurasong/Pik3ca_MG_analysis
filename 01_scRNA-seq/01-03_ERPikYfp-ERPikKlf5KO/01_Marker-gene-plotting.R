#####################################################################################
# UMAP feature plots of epithelial marker genes (ER_Pik vs ER_Pik Klf5 KO)
#
# Description:
#   Generates UMAP-based expression overlays for selected epithelial marker genes
#   in two annotated epithelial Seurat objects: ER Pik Yfp (control oncogene) and
#   ER Pik Klf5 knockout. The script extracts UMAP embeddings and log-normalized
#   RNA expression values, then creates ggplot-based feature maps using a
#   consistent viridis colour gradient and minimalist theme for figure assembly.
#   These panels correspond to Fig. 5e (ER Pik Yfp) and Fig. 5f (ER Pik Klf5 KO).
#
# Workflow overview:
#     1. Load annotated epithelial Seurat objects (ER_Pik and ER_Pik_Klf5KO)
#     2. Define a shared continuous colour palette (viridis; reversed)
#     3. For each object:
#         a. Extract UMAP coordinates
#         b. Merge UMAP coordinates with RNA assay expression matrix (`data`)
#         c. Plot gene expression on UMAP using custom ggplot FeaturePlot-like function
#     4. Export/arrange plots downstream for figure panels
#
# Inputs:
#   - ERPik_Epith_annotate.rds             : ER Pik Yfp epithelial Seurat object
#   - ERPik_Klf5_Epith_annotate.rds        : ER Pik Klf5 KO epithelial Seurat object
#
# Outputs:
#   - UMAP expression plots for marker genes used in:
#       Fig. 5e (ER Pik Yfp) and Fig. 5f (ER Pik Klf5 KO)
#   - A list of ggplot objects (one per gene, per condition) for downstream assembly
#
# Marker genes plotted:
#   - Krt8   : luminal epithelial marker
#   - Esr1   : ER+ luminal marker
#   - Kit    : progenitor-associated marker
#   - Krt14  : basal/myoepithelial marker
#
# Dependencies:
#   Seurat, tidyverse, dplyr, data.table, ggplot2, viridis,
#   RColorBrewer, plotly, grid
#
# Notes:
#   - Expression is taken from `object[["RNA"]]$data` (log-normalized values).
#   - Legends and axis elements are removed to facilitate multi-panel figure layout.
#   - Ensure gene symbols exist in each object before plotting (some markers may be absent).
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

erpik <- readRDS("ERPik_Epith_annotate.rds")
klf5 <- readRDS("ERPik_Klf5_Epith_annotate.rds")

# Colour palette

myPalette <- viridis(n = 10, option = "C", direction = -1)
sc <- scale_colour_gradientn(colours = rev(myPalette))

# Plotting of ER Pik Yfp (Fig. 5e)

## Set-up

umap_coord <- as.data.frame(erpik@reductions$umap@cell.embeddings)
erpik_gene <- merge.data.frame(umap_coord,t(erpik[["RNA"]]$data),by=0)
row.names(erpik_gene) <- erpik_gene$Row.names

## Code for visualisation

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- erpik_gene[,regulon]
  regulon_plot <- ggplot(data=erpik_gene) + geom_point(aes(x=umap_1,y=umap_2,colour=regulon_d), size=1)
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

gene_list <- c("Krt8", "Esr1", "Kit", "Krt14")
lapply(gene_list, erpik)

# Plotting of ER Pik Klf5 KO (Fig. 5f)

## Set-up

umap_coord <- as.data.frame(klf5@reductions$umap@cell.embeddings)
klf5_gene <- merge.data.frame(umap_coord,t(klf5[["RNA"]]$data),by=0)
row.names(klf5_gene) <- klf5_gene$Row.names

## Code for visualisation

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- klf5_gene[,regulon]
  regulon_plot <- ggplot(data=klf5_gene) + geom_point(aes(x=umap_1,y=umap_2,colour=regulon_d), size=1)
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

gene_list <- c("Krt8", "Esr1", "Kit", "Krt14")
lapply(gene_list, klf5)
