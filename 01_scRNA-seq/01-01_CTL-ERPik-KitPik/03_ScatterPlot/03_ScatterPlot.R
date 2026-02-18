#####################################################################################
# Gene-set detection ratios and adjusted lineage proportions (Fig. 2i–j, 2m; Ext. Fig. 4d)
#
# Description:
#   Quantifies per-cell detection ratios for predefined marker gene sets
#   (BC/Myo, ER+ LCs, ER- LCs) using raw RNA counts from the integrated epithelial
#   Seurat object. To correct for sequencing depth/complexity effects, detection
#   fractions are adjusted using robust linear models (rlm) against nFeature_RNA.
#   The adjusted proportions are then:
#     (1) visualized as scatter plots comparing lineage programs in HY states
#         (Early HY: Fig. 2i; Late HY: Fig. 2j),
#     (2) used to compare ER+ vs ER- LC programs (Fig. 2m),
#     (3) projected back onto the UMAP embedding as split FeaturePlots
#         (Extended Data Fig. 4d).
#
# Workflow overview:
#     1. Load integrated epithelial Seurat object
#     2. Load marker gene lists (BC/Myo, ER+ LC, ER- LC)
#     3. Extract RNA count matrix and selected metadata (incl. nFeature_RNA, sample, cell_type)
#     4. Compute per-cell detection fractions:
#          fraction = (# marker genes detected >0) / (size of marker set)
#     5. Adjust detection fractions for nFeature_RNA using robust linear regression (rlm)
#          adjusted_prop = mean(fraction) + residuals(rlm(fraction ~ nFeature_RNA))
#        performed for:
#          a) BC vs ER+ LC fractions
#          b) ER+ LC vs ER- LC fractions
#     6. Generate scatter plots of adjusted proportions for specific subsets:
#          - CTL reference cells + Early HY (ER_Pik)  → Fig. 2i
#          - CTL reference cells + Late HY (ER_Pik)   → Fig. 2j
#          - CTL ER- LC subsets                       → Fig. 2m
#     7. Replace/attach adjusted proportions to meta.data and generate split UMAP panels
#        using a custom wrapper (auc_panel) with shared color scaling and PDF export.
#
# Inputs:
#   - Epithelial_Integrated_Annot.rds
#       Integrated epithelial Seurat object with UMAP and annotations (cell_type, sample)
#   - marker_genes/BC_markers.rds
#   - marker_genes/ER+_LC_markers.rds
#   - marker_genes/ER-_LC_markers.rds
#
# Outputs:
#   - Scatter plots (adjusted lineage proportions):
#       Fig. 2i : Early HY vs CTL reference (BC/Myo vs ER+ LC)
#       Fig. 2j : Late HY  vs CTL reference (BC/Myo vs ER+ LC)
#       Fig. 2m : ER+ LC vs ER- LC (CTL subsets)
#   - Split UMAP panels (PDF export; Extended Data Fig. 4d):
#       adjusted_bc_prop_AUC_panel.pdf
#       adjusted_lcp_prop_AUC_panel.pdf
#       adjusted_lcm_prop_AUC_panel.pdf
#
# Dependencies:
#   Seurat, SingleCellExperiment, ggplot2, Matrix, dplyr, patchwork,
#   clustree, tidyverse, scales, MASS
#
# Notes:
#   - Fractions are based on detection (non-zero counts), not expression magnitude.
#   - Robust regression reduces the influence of outliers when correcting for
#     sequencing depth/complexity (nFeature_RNA).
#   - This script overwrites `seuset_adj@meta.data` with merged_adjust; ensure
#     rownames match cell barcodes exactly to avoid misalignment.
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
seuset <- readRDS("Epithelial_Integrated_Annot.rds") 

# Marker genes
bc <- readRDS("marker_genes/BC_markers.rds")
lcp <- readRDS("marker_genes/ER+_LC_markers.rds")
lcm <- readRDS("marker_genes/ER-_LC_markers.rds")

# Get information from meta.data
meta <- seuset@meta.data[,c(3, 6, 8, 11)]

# Calculate ratio
c_seuset <- as.matrix(seuset[["RNA"]]$counts) 

bc_seuset <- c_seuset[rownames(c_seuset) %in% bc, ]
summary(colSums(bc_seuset != 0))
fraction_bc <- colSums(bc_seuset != 0)/nrow(bc_seuset)
bc_seuset_up <- rbind(bc_seuset, fraction_bc)

lcp_seuset <- c_seuset[rownames(c_seuset) %in% lcp, ]
summary(colSums(lcp_seuset != 0))
fraction_lcp <- colSums(lcp_seuset != 0)/nrow(lcp_seuset)
lcp_seuset_up <- rbind(lcp_seuset, fraction_lcp)

lcm_seuset <- c_seuset[rownames(c_seuset) %in% lcm, ]
summary(colSums(lcm_seuset != 0))
fraction_lcm <- colSums(lcm_seuset != 0)/nrow(lcm_seuset)
lcm_seuset_up <- rbind(lcm_seuset, fraction_lcm)

# Ratio adjustment
## BC and ER+ LCs

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

bc_l <- rlm(fraction_bc ~ nFeature_RNA, data=bc_lcp_mer)
lcp_l <- rlm(fraction_lcp ~ nFeature_RNA, data=bc_lcp_mer)

bc_fitted <- mean(bc_lcp_mer$fraction_bc) + bc_l$residuals[rownames(bc_lcp_mer)]
lcp_fitted <- mean(bc_lcp_mer$fraction_lcp) + lcp_l$residuals[rownames(bc_lcp_mer)]

bc_lcp_mer[rownames(bc_lcp_mer), "adjusted_bc_prop"] <- bc_fitted
bc_lcp_mer[rownames(bc_lcp_mer), "adjusted_lcp_prop"] <- lcp_fitted

## LC ER+ - LC ER- 

lcm_frac <- as.matrix(lcm_seuset_up[nrow(lcm_seuset_up),])
colnames(lcm_frac) <- "fraction_lcm"

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

# Plotting scatter plot
## BC - ER+ LCs in Early HY (Fig. 2i)

ctl <- bc_lcp_mer[bc_lcp_mer$cell_type %in% c("Myo", 'ER+ LCs Sca1', 'ER+ LCs Foxa1'),]
ctl <- ctl[ctl$sample == "CTL", ]

early_hy <- bc_lcp_mer[bc_lcp_mer$cell_type == "Early HY",]
early_hy <- early_hy[early_hy$sample == "ER_Pik", ]

bc_sel <- rbind(ctl, early_hy)

ggplot(data=bc_sel, aes(x=adjusted_bc_prop, y=adjusted_lcp_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#BDB76B", "#800000", "#F08080", "#1E90FF")) +
  xlab("Adjusted Myo proportion") +
  ylab("Adjusted LC ER+ proportion") +
  geom_abline(slope = 1) +
  ylim(0, 1) +
  xlim(0, 1) + 
  NoLegend() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## BC - ER+ LCs in Late HY (Fig. 2j)

late_hy <- bc_lcp_mer[bc_lcp_mer$cell_type == "Late HY",]
late_hy <- late_hy[late_hy$sample == "ER_Pik", ]

late_sel <- rbind(ctl, late_hy)

ggplot(data=late_sel, aes(x=adjusted_bc_prop, y=adjusted_lcp_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#800000", "#F08080", "#32CD32", "#1E90FF")) +
  xlab("Adjusted Myo proportion") +
  ylab("Adjusted LC ER+ proportion") +
  geom_abline(slope = 1) +
  ylim(0, 1) +
  xlim(0, 1) + 
  NoLegend() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

## ER+ LCs - ER- LCs (Fig. 2m)

bc_sel <- lcp_lcm_mer[lcp_lcm_mer$cell_type %in% c('ER- LCs Itgb3low', "ER- LCs Itgb3high"),]
bc_sel <- bc_sel[bc_sel$sample == "CTL", ]

ggplot(data=bc_sel, aes(x=adjusted_lcp_prop, y=adjusted_lcm_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#FFA500", "#2F4F4F")) +
  xlab("Adjusted LC ER+ proportion") +
  ylab("Adjusted LC ER- proportion") +
  geom_abline(slope = 1) +
  ylim(0, 1) +
  xlim(0, 1) + 
  NoLegend() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

# Plotting as split UMAP
## Embedding

bc_lcp_adjust <- bc_lcp_mer[, 7:8]
bc_lcm_adjust <- bc_lcm_mer[, c(5:8)]
merge_adjust <- merge(bc_lcp_adjust, bc_lcm_adjust, by=0)
rownames(merge_adjust) <- merge_adjust$Row.names
merge_adjust <- merge_adjust[, -1]

merge_adjust <- merge_adjust[, c(1, 2, 6, 3, 4)]
names(merge_adjust) <- c("adjusted_bc_prop", "adjusted_lcp_prop", "adjusted_lcm_prop", "cell_type", "sample")

seuset_adj <- seuset
seuset_adj@meta.data <- merge_adjust

## Plotting 

auc_panel <- function(obj, feature,
                      split.by = "sample",
                      outdir = "AUC_plot_split",
                      palette_option = "C",
                      lims = NULL,
                      reduction = "umap") {
  # feature must exist in obj@meta.data
  if (!feature %in% colnames(obj@meta.data)) {
    stop(sprintf("`%s` not found in obj@meta.data", feature))
  }
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # limits for consistent legend
  if (is.null(lims)) {
    vals <- obj@meta.data[[feature]]
    lims <- range(vals, na.rm = TRUE, finite = TRUE)
  }
  
  pal <- viridis(10, option = palette_option, direction = -1)
  sc  <- scale_colour_gradientn(
    colours = rev(pal),
    limits  = lims,
    oob     = squish,
    na.value = "grey90"
  )
  
  # split FeaturePlot
  fp_list <- FeaturePlot(
    obj,
    features   = feature,
    split.by   = split.by,
    keep.scale = "all",
    order      = TRUE,
    reduction  = reduction,
    combine    = FALSE
  )
  fp_list <- lapply(fp_list, `+`, sc)
  
  # get levels of split.by for proper labels
  lvls <- levels(factor(obj@meta.data[[split.by]]))
  
  # add titles and remove UMAP axes
  fp_list <- Map(function(p, lab) {
    p +
      ggtitle(lab) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank(),
        panel.grid = element_blank()
      )
  }, fp_list, lvls)
  
  # combine panels with shared legend and big title
  p <- wrap_plots(fp_list) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  p <- p + plot_annotation(
    title = feature,
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )
  
  # save
  ggsave(
    filename = file.path(outdir, paste0(feature, "_AUC_panel.pdf")),
    plot     = p, height = 4, width = 10, limitsize = FALSE
  )
  
  return(p)
}

## Extended Data Fig. 4d
auc_panel(seuset_adj, "adjusted_bc_prop")
auc_panel(seuset_adj, "adjusted_lcp_prop")
auc_panel(seuset_adj, "adjusted_lcm_prop")

