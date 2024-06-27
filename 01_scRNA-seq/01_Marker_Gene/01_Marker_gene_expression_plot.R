#####################################################################################
# This code is applied on K8Pik data.
# All parameters for expression plotting are exactly same for CTL and ERPik data.
# You just need to change the object for the individual data.
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

seuset <- readRDS("K8Pik_annotated.rds") # You can change here to plot the other datasets.

umap_coord <- as.data.frame(seuset@reductions$umap@cell.embeddings)

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
  
  ggsave(paste(x,"_GeneExp_K8Pik.pdf",sep=""), width = 10, height = 10, units = "cm")
}

gene_list <- read.delim("gene_list.txt", h=F)
gene_list <- gene_list$V1

lapply(gene_list, umap_color_scaled)