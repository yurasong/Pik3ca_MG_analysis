#####################################################################################
# Differential expression analysis and volcano plots for HY state transitions
#
# Description:
#   Performs differential gene expression (DEG) analysis on the integrated
#   epithelial Seurat object to identify genes upregulated during transition
#   into hybrid (HY) states. Specifically:
#     (1) Early HY vs ER+ luminal cells (ER+ LCs Sca1 + ER+ LCs Foxa1)
#     (2) Late HY vs Early HY
#   Results are visualised as volcano plots highlighting significantly
#   upregulated genes (adjusted P-value < 0.01). These panels correspond to
#   Fig. 4a (Early HY vs ER+ LCs) and Fig. 4b (Late HY vs Early HY).
#
# Workflow overview:
#     1. Load integrated epithelial Seurat object and set DefaultAssay to RNA
#     2. Identify upregulated markers using FindMarkers (only.pos = TRUE)
#     3. Add significance labels based on adjusted p-value threshold
#     4. Compute -log10(adjusted p-value) for volcano plot y-axis
#     5. Generate volcano plots for each comparison
#
# Inputs:
#   - Epithelial_Integrated_Annot.rds
#       Integrated epithelial Seurat object with `cell_type` identities
#
# Outputs:
#   - Volcano plot: Early HY vs ER+ LCs (Fig. 4a)
#   - Volcano plot: Late HY vs Early HY (Fig. 4b)
#   - (Optional; recommended) export DEG tables for downstream plotting/reporting
#
#
# Dependencies:
#   Seurat, ggplot2, tidyverse, patchwork, data.table, magrittr
#
# Notes:
#   - This script currently plots all points; labeling top genes can be added
#     downstream (e.g., ggrepel).
#   - Ensure the identities (Idents) are set to the intended `cell_type` field
#     prior to FindMarkers if not already done in the saved object.
#####################################################################################

# Library
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

# Data

seuset <- readRDS("Epithelial_Integrated_Annot.rds")
DefaultAssay(seuset) <- "RNA"

# DEG analysis
## GeneUp Early HY vs ER+ LC

up_early <- FindMarkers(seuset, 
                        ident.1="Early HY",
                        ident.2=c("ER+ LCs Sca1", "ER+ LCs Foxa1"),
                        only.pos=T,
                        min.pct=0.25, logfc.threshold=0.5)

## GeneUp Late HY vs Early HY

late_early <- FindMarkers(seuset, 
                        ident.1="Late HY",
                        ident.2="Early HY",
                        only.pos=T,
                        min.pct=0.25, logfc.threshold=0.5)

# Visualisation
## Volcano Plot of GeneUp Early HY vs ER+ LC (Fig. 4a)

up_early$sig <- "no"
up_early$sig[up_early$p_val_adj < 0.01] <- "Significant"
up_early$sig[up_early$p_val_adj >= 0.01] <- "Non_Significant"

up_early$logP <- -log10(up_early$p_val_adj)

ggplot(data=earlyhy_vs_erp, aes(x=avg_log2FC, y=logP, col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal() +
  ylab('-log(adjusted_Pvalue)')

## Volcano Plot of GeneUp Late HY vs Early HY (Fig. 4b)

late_early$sig <- "no"
late_early$sig[late_early$p_val_adj < 0.01] <- "Significant"
late_early$sig[late_early$p_val_adj >= 0.01] <- "Non_Significant"

up_elate_earlyarly$logP <- -log10(late_early$p_val_adj)

ggplot(data=latehy_vs_earlyhy, aes(x=avg_log2FC, y=logP, col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal() +
  ylab('-log(adjusted_Pvalue)')
