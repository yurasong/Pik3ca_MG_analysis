#######################################################################################
# ScatterPlot_KitPik.R
# This script will generate the scatter plot of BC-LC marker genes on Kit Pik Yfp data.
# This code will generate plots in Fig. 2i.
#######################################################################################

# Library
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(dplyr)
library(patchwork)
library(clustree)
library(tidyverse)
library(scales)
library(MASS)

# Data
seuset <- readRDS("KitPik_annotated.rds")

# Get information from meta.data
meta <- seuset@meta.data[,c(3, 9)]

## Marker genes
bc <- readRDS("marker_genes/BC_markers.rds")
lcm <- readRDS("marker_genes/ER-_LC_markers.rds")

## Get the count table
c_seuset <- as.matrix(seuset@assays[["RNA"]]$counts) # Seurat assay V5

### working on BC
bc_seuset <- c_seuset[rownames(c_seuset) %in% bc, ]
fraction_bc <- colSums(bc_seuset != 0)/nrow(bc_seuset)
bc_seuset_up <- rbind(bc_seuset, fraction_bc)

### working on LC ER-
lcm_seuset <- c_seuset[rownames(c_seuset) %in% lcm, ]
fraction_lcm <- colSums(lcm_seuset != 0)/nrow(lcm_seuset)
lcm_seuset_up <- rbind(lcm_seuset, fraction_lcm)

# Plotting

## BC and LC ER-

bc_frac <- as.matrix(bc_seuset_up[nrow(bc_seuset_up),])
colnames(bc_frac) <- "fraction_bc"
head(bc_frac)

lcm_frac <- as.matrix(lcm_seuset_up[nrow(lcm_seuset_up),])
colnames(lcm_frac) <- "fraction_lcm"
head(lcm_frac)

bc_lcm <- merge(bc_frac, lcm_frac, by="row.names")
rownames(bc_lcm) <- bc_lcm$Row.names
bc_lcm <- bc_lcm[, 2:3]
head(bc_lcm)

bc_lcm_mer <- merge(bc_lcm, meta, by="row.names")

rownames(bc_lcm_mer) <- bc_lcm_mer$Row.names
bc_lcm_mer <- bc_lcm_mer[, 2:ncol(bc_lcm_mer)]
head(bc_lcm_mer)

bc_l <- rlm(fraction_bc ~ nFeature_RNA, data=bc_lcm_mer)
lcm_l <- rlm(fraction_lcm ~ nFeature_RNA, data=bc_lcm_mer)

bc_fitted <- mean(bc_lcm_mer$fraction_bc) + bc_l$residuals[rownames(bc_lcm_mer)]
lcm_fitted <- mean(bc_lcm_mer$fraction_lcm) + lcm_l$residuals[rownames(bc_lcm_mer)]

bc_lcm_mer[rownames(bc_lcm_mer), "adjusted_bc_prop"] <- bc_fitted
bc_lcm_mer[rownames(bc_lcm_mer), "adjusted_lcm_prop"] <- lcm_fitted

head(bc_lcm_mer)

### Plotting
bc_sel <- bc_lcm_mer[bc_lcm_mer$cell_type %in% c("Myo", "ER- LCs"),]

ggplot(data=bc_sel, aes(x=adjusted_bc_prop, y=adjusted_lcm_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#2E8B57", "#1E90FF")) +
  xlab("Adj BC prop") +
  ylab("Adj ER- prop") +
  geom_abline(slope = 1) +
  ylim(0.1, 1) +
  xlim(0.1, 1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none", 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )