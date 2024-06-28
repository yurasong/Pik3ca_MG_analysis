#####################################################################################
# This code is applied on CTL data.
#####################################################################################

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

seuset <- readRDS("Control_annotated.rds")

# Get information from meta.data
meta <- seuset@meta.data[,c(3, 8)]

## Marker genes
bc <- readRDS("marker_genes/BC_markers.rds")
lcp <- readRDS("marker_genes/ER+_LC_markers.rds")
lcm <- readRDS("marker_genes/ER-_LC_markers.rds")

## Get the count table
c_seuset <- as.matrix(seuset@assays[["RNA"]]@counts)

### working on BC
bc_seuset <- c_seuset[rownames(c_seuset) %in% bc, ]
fraction_bc <- colSums(bc_seuset != 0)/nrow(bc_seuset)
bc_seuset_up <- rbind(bc_seuset, fraction_bc)

### working on LC ER+
lcp_seuset <- c_seuset[rownames(c_seuset) %in% lcp, ]
fraction_lcp <- colSums(lcp_seuset != 0)/nrow(lcp_seuset)
lcp_seuset_up <- rbind(lcp_seuset, fraction_lcp)

### working on LC ER-
lcm_seuset <- c_seuset[rownames(c_seuset) %in% lcm, ]
fraction_lcm <- colSums(lcm_seuset != 0)/nrow(lcm_seuset)
lcm_seuset_up <- rbind(lcm_seuset, fraction_lcm)

# Plotting

## BC and LC ER+ 
bc_frac <- as.matrix(bc_seuset_up[nrow(bc_seuset_up),])
colnames(bc_frac) <- "fraction_bc"

lcp_frac <- as.matrix(lcp_seuset_up[nrow(lcp_seuset_up),])
colnames(lcp_frac) <- "fraction_lcp"

bc_lcp <- merge(bc_frac, lcp_frac, by="row.names")
rownames(bc_lcp) <- bc_lcp$Row.names
bc_lcp <- bc_lcp[, 2:3]

bc_lcp_mer <- merge(bc_lcp, meta, by="row.names")

rownames(bc_lcp_mer) <- bc_lcp_mer$Row.names
bc_lcp_mer <- bc_lcp_mer[, 2:ncol(bc_lcp_mer)]
head(bc_lcp_mer)

bc_l <- rlm(fraction_bc ~ nFeature_RNA, data=bc_lcp_mer)
lcp_l <- rlm(fraction_lcp ~ nFeature_RNA, data=bc_lcp_mer)

bc_fitted <- mean(bc_lcp_mer$fraction_bc) + bc_l$residuals[rownames(bc_lcp_mer)]
lcp_fitted <- mean(bc_lcp_mer$fraction_lcp) + lcp_l$residuals[rownames(bc_lcp_mer)]
bc_lcp_mer[rownames(bc_lcp_mer), "adjusted_bc_prop"] <- bc_fitted
bc_lcp_mer[rownames(bc_lcp_mer), "adjusted_lcp_prop"] <- lcp_fitted

### Plotting
bc_sel <- bc_lcp_mer[bc_lcp_mer$cell_type %in% c("BC", "LC_ER+"),]
ggplot(data=bc_sel, aes(x=adjusted_bc_prop, y=adjusted_lcp_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#FF6347", "#1E90FF")) +
  xlab("Adj BC prop") +
  ylab("Adj ER+ prop") +
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

## BC and LC ER-
lcm_frac <- as.matrix(lcm_seuset_up[nrow(lcm_seuset_up),])
colnames(lcm_frac) <- "fraction_lcm"

bc_lcm <- merge(bc_frac, lcm_frac, by="row.names")
rownames(bc_lcm) <- bc_lcm$Row.names
bc_lcm <- bc_lcm[, 2:3]

bc_lcm_mer <- merge(bc_lcm, meta, by="row.names")
rownames(bc_lcm_mer) <- bc_lcm_mer$Row.names
bc_lcm_mer <- bc_lcm_mer[, 2:ncol(bc_lcm_mer)]

bc_l <- rlm(fraction_bc ~ nFeature_RNA, data=bc_lcm_mer)
lcm_l <- rlm(fraction_lcm ~ nFeature_RNA, data=bc_lcm_mer)

bc_fitted <- mean(bc_lcm_mer$fraction_bc) + bc_l$residuals[rownames(bc_lcm_mer)]
lcm_fitted <- mean(bc_lcm_mer$fraction_lcm) + lcm_l$residuals[rownames(bc_lcm_mer)]
bc_lcm_mer[rownames(bc_lcm_mer), "adjusted_bc_prop"] <- bc_fitted
bc_lcm_mer[rownames(bc_lcm_mer), "adjusted_lcm_prop"] <- lcm_fitted

### Plotting
bc_sel <- bc_lcm_mer[bc_lcm_mer$cell_type %in% c("BC", "LC_ER-_1", "LC_ER-_2"),]

ggplot(data=bc_sel, aes(x=adjusted_bc_prop, y=adjusted_lcm_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#2E8B57", "#1E90FF", "#FFA500")) +
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

## LC ER+ and LC ER-

lcp_lcm <- merge(lcp_frac, lcm_frac, by="row.names")
rownames(lcp_lcm) <- lcp_lcm$Row.names
lcp_lcm <- lcp_lcm[, 2:3]

lcp_lcm_mer <- merge(lcp_lcm, meta, by="row.names")
rownames(lcp_lcm_mer) <- lcp_lcm_mer$Row.names
lcp_lcm_mer <- lcp_lcm_mer[, 2:ncol(lcp_lcm_mer)]

lcp_l <- rlm(fraction_lcp ~ nFeature_RNA, data=lcp_lcm_mer)
lcm_l <- rlm(fraction_lcm ~ nFeature_RNA, data=lcp_lcm_mer)

lcp_fitted <- mean(lcp_lcm_mer$fraction_lcp) + lcp_l$residuals[rownames(lcp_lcm_mer)]
lcm_fitted <- mean(lcp_lcm_mer$fraction_lcm) + lcm_l$residuals[rownames(lcp_lcm_mer)]
lcp_lcm_mer[rownames(lcp_lcm_mer), "adjusted_lcp_prop"] <- lcp_fitted
lcp_lcm_mer[rownames(lcp_lcm_mer), "adjusted_lcm_prop"] <- lcm_fitted

### Plotting
bc_sel <- lcp_lcm_mer[lcp_lcm_mer$cell_type %in% c("LC_ER+", "LC_ER-_1", "LC_ER-_2"),]

ggplot(data=bc_sel, aes(x=adjusted_lcp_prop, y=adjusted_lcm_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#FF6347", "#2E8B57", "#FFA500")) +
  xlab("Adj ER+ prop") +
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

