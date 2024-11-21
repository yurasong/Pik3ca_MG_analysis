#####################################################################################
# This code is applied on Integrated data.
#####################################################################################

# Library
library(Signac)
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
seuset <- readRDS("Annot_multiome_integrated.rds")

# Get information from meta.data
meta <- seuset@meta.data[,c(3, 17)]
names(meta) <- c("nFeature_RNA", "cell_type")
head(meta)

# Marker genes
bc <- readRDS("05_marker_genes/BC_markers.rds")
lcp <- readRDS("05_marker_genes/ER+_LC_markers.rds")
lcm <- readRDS("05_marker_genes/ER-_LC_markers.rds")

# Calculate the fraction and add on the dataset

## working on BC
bc_seuset <- c_seuset[rownames(c_seuset) %in% bc, ]
fraction_bc <- colSums(bc_seuset != 0)/nrow(bc_seuset)
bc_seuset_up <- rbind(bc_seuset, fraction_bc)

## working on LC ER+
lcp_seuset <- c_seuset[rownames(c_seuset) %in% lcp, ]
fraction_lcp <- colSums(lcp_seuset != 0)/nrow(lcp_seuset)
lcp_seuset_up <- rbind(lcp_seuset, fraction_lcp)

## working on LC ER-
lcm_seuset <- c_seuset[rownames(c_seuset) %in% lcm, ]
fraction_lcm <- colSums(lcm_seuset != 0)/nrow(lcm_seuset)
lcm_seuset_up <- rbind(lcm_seuset, fraction_lcm)

# Plotting

## BC and LC ER+ 
bc_frac <- as.matrix(bc_seuset_up[nrow(bc_seuset_up),])
colnames(bc_frac) <- "fraction_bc"
head(bc_frac)

lcp_frac <- as.matrix(lcp_seuset_up[nrow(lcp_seuset_up),])
colnames(lcp_frac) <- "fraction_lcp"
head(lcp_frac)

bc_lcp <- merge(bc_frac, lcp_frac, by="row.names")
rownames(bc_lcp) <- bc_lcp$Row.names
bc_lcp <- bc_lcp[, 2:3]
head(bc_lcp)

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

### Reference
bc_lcp_sel <- bc_lcp_mer[bc_lcp_mer$cell_type %in% c("LC_ER+_Klf5_specific", "LC_ER+", "Myoepith"), ]

ggplot(data=bc_lcp_sel, aes(x=adjusted_bc_prop, y=adjusted_lcp_prop, color = cell_type)) +
  geom_point() +
  xlab("Adjusted BC propertion") +
  ylab("Adjusted LC ER+ propertion") +
  geom_abline(slope = 1) +
  ylim(0, 0.8) +
  xlim(0, 0.8) +
  scale_colour_manual(values=c("#00A9FF", "#C77CFF", "#800080")) +
  theme(legend.position = "none")

### BCs
bc_lcp_sel <- bc_lcp_mer[bc_lcp_mer$cell_type %in% c("Immature_BC"), ]

ggplot(data=bc_lcp_sel, aes(x=adjusted_bc_prop, y=adjusted_lcp_prop, color = cell_type)) +
  geom_point() +
  xlab("Adjusted BC propertion") +
  ylab("Adjusted LC ER+ propertion") +
  geom_abline(slope = 1) +
  ylim(0, 0.8) +
  xlim(0, 0.8) +
  scale_colour_manual(values=c("#00BE67")) +
  theme(legend.position = "none")

### HY ER+/BC
bc_lcp_sel <- bc_lcp_mer[bc_lcp_mer$cell_type %in% c("HY_ER+/BC"), ]

ggplot(data=bc_lcp_sel, aes(x=adjusted_bc_prop, y=adjusted_lcp_prop, color = cell_type)) +
  geom_point() +
  xlab("Adjusted BC propertion") +
  ylab("Adjusted LC ER+ propertion") +
  geom_abline(slope = 1) +
  ylim(0, 0.8) +
  xlim(0, 0.8) +
  scale_colour_manual(values=c("#CD9600")) +
  theme(legend.position = "none")

## BC and LC ER-
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

### Reference
bc_lcm_sel <- bc_lcm_mer[bc_lcm_mer$cell_type %in% c("LC_ER-", "Myoepith"), ]

ggplot(data=bc_lcm_sel, aes(x=adjusted_bc_prop, y=adjusted_lcm_prop, color = cell_type)) +
  geom_point() +
  xlab("Adjusted BC propertion") +
  ylab("Adjusted LC ER- propertion") +
  geom_abline(slope = 1) +
  ylim(0, 0.8) +
  xlim(0, 0.8) +
  scale_colour_manual(values=c("#00BFC4", "#C77CFF")) +
  theme(legend.position = "none")

### HY BC/ER-
bc_lcm_sel <- bc_lcm_mer[bc_lcm_mer$cell_type %in% c("HY BC/ER-"), ]

ggplot(data=bc_lcm_sel, aes(x=adjusted_bc_prop, y=adjusted_lcm_prop, color = cell_type)) +
  geom_point() +
  xlab("Adjusted BC propertion") +
  ylab("Adjusted LC ER- propertion") +
  geom_abline(slope = 1) +
  ylim(0, 0.8) +
  xlim(0, 0.8) +
  scale_colour_manual(values=c("#F8766D")) +
  theme(legend.position = "none")

