#####################################################################################
# 01_Marker_gene_expression_plot.R
# Change the `seuset <- readRDS(...)` line to point at the desired RDS.
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

# Load data
seuset <- readRDS("seurat_object.rds") # swap filename here for other datasets

# extract UMAP coords
umap_coord <- as.data.frame(seuset@reductions$umap@cell.embeddings)

# merge expression with coords
seuset_gene <- merge.data.frame(umap_coord,t(as.data.frame(seuset@assays$RNA@data)),by=0)
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
}

# --- Run for each marker ----------------------------------------------------
gene_list <- read.delim("gene_list.txt", header = FALSE)$V1

lapply(gene_list, umap_color_scaled)