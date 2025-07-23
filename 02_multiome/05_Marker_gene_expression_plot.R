#####################################################################################
# File: 05_Marker_gene_expression_plot.R
#
# Description:
#   Generates UMAP plots colored by expression of selected marker genes
#   on the integrated multiome Seurat object.
#
# Inputs:
#   - Multiome_integrated_annot.rds : integrated Seurat object with UMAP reductions
#
# Outputs: Plots on Fig. S5a
#   - PDF files named "<gene>_GeneExp_rescale-Integrated.pdf" for each gene in `gene_list`
#
# Dependencies:
#   Seurat, Signac, RColorBrewer, dplyr, plotly, grid, data.table, tidyverse, viridis
#
# Notes:
#   - Expects the “RNA” assay to contain normalized expression data.
#   - UMAP coordinates are taken from `seuset@reductions$umap@cell.embeddings`.
#####################################################################################


# Library
library(Seurat)
library(Signac)
library(RColorBrewer)
library(dplyr)
library(plotly)
library(grid)
library(data.table)
library(tidyverse)
library(viridis)

# Data

seuset <- readRDS("Multiome_integrated_annot.rds")

umap_coord <- as.data.frame(seuset@reductions$umap@cell.embeddings)

seuset_gene <- merge.data.frame(umap_coord,t(as.data.frame(seuset[[RNA]]$data)),by=0) # V5 assay
row.names(seuset_gene) <- seuset_gene$Row.names

# Plotting

myPalette <- viridis(n = 10, option = "C", direction = -1)
sc <- scale_colour_gradientn(colours = rev(myPalette))

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
  
  ggsave(paste(x,"_GeneExp_rescale-Integrated.pdf",sep=""), width = 10, height = 10, units = "cm")
}

gene_list <- c("Krt5", "Krt14", "Acta2", "Myh11", "Mylk",
               "Krt8", "Foxa1", "Prlr", "Esr1", "Prom1",
               "Elf5", "Kit", "Cd14", "Itga2", "Csn3",
               "Il33", "Krt6a", "Klf5") # list of genes we checked for annotation 

lapply(gene_list, umap_color_scaled)