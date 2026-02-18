#######################################################################################
# 06_embed_data_on_seurat.R
# This script will embed average AUC score into the Seurat object as a new assay.
# Plotting of AUC value will generate the graphs on Fig. 5d-g.
#######################################################################################

# Preparations

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(dplyr)
library(patchwork)
library(clustree)
library(tidyverse)
library(pheatmap)

# Prepare the datasets

## Average result

scenic_results <- read.table("AUC_Average.csv", header=T, sep=",", row.names=1) # Avarage AUC value from 04_merge_results.py
colnames(scenic_results) <- sub("...$","",colnames(scenic_results))

## Seurat proceeded object

seuset <- readRDS("seurat_object.rds")

# Embedding on Seurat

# AUC matrix as the assay of seurat object
## To put the AUC matrix as additional assay on seurat object

seuset_auc <- seuset # To avoid issues on over-writing

auc_mtx <- t(scenic_results)
aucs <- CreateAssayObject(counts = auc_mtx)
seuset_auc[["AUC"]] <- aucs

DefaultAssay(seuset_auc) <- "AUC"

# Plotting AUC value

umap_coord <- as.data.frame(seuset_auc@reductions$umap@cell.embeddings)
seuset_auc_matrix <- merge.data.frame(umap_coord,t(as.data.frame(seuset_auc@assays$AUC@data)),by=0)
row.names(seuset_auc_matrix) <- seuset_auc_matrix$Row.names
seuset_auc_matrix <- seuset_auc_matrix[, -1]

# Plotting

myPalette <- viridis(n = 10, option = "C", direction = -1)
sc <- scale_colour_gradientn(colours = rev(myPalette))

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- seuset_auc_matrix[,regulon]
  regulon_plot <- ggplot(data=seuset_auc_matrix) + geom_point(aes(x=umap_1,y=umap_2,colour=regulon_d), size=1)
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

umap_color_scaled("Klf5")