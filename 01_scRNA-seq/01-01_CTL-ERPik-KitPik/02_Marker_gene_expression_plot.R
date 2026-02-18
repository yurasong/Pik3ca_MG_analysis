#####################################################################################
# FeaturePlot visualisation of key epithelial marker genes across samples (Fig. 2c–h)
#
# Description:
#   Generates split-by-sample FeaturePlot panels for canonical epithelial marker
#   genes using the integrated epithelial Seurat object. Expression is displayed
#   on the UMAP embedding, split by experimental condition (`sample`), and a
#   consistent viridis colour gradient is applied for publication-ready figure
#   panels. These plots correspond to Fig. 2c–h.
#
# Workflow overview:
#     1. Load integrated epithelial Seurat object
#     2. Define a shared continuous viridis colour palette (reversed)
#     3. Generate FeaturePlot panels split by `sample`
#     4. Apply consistent colour scaling (and quantile cutoffs where indicated)
#
# Inputs:
#   - Epithelial_Integrated_Annot.rds
#       Integrated epithelial Seurat object containing UMAP embedding and
#       metadata column `sample`
#
# Dependencies:
#   Seurat, tidyverse, dplyr, data.table, viridis,
#   RColorBrewer, plotly, grid
#
# Notes:
#   - Quantile-based cutoffs are used for selected genes to reduce the influence
#     of outliers and improve contrast.
#   - Ensure DefaultAssay is set appropriately (typically "RNA") before plotting.
#   - For consistency, assign FeaturePlot output to `p <- ...` before applying
#     patchwork scaling (e.g., `p & scale_colour_gradientn(...)`).
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
seuset <- readRDS("Epithelial_Integrated_Annot.rds") 

# Palette
myPalette <- viridis(n = 10, option = "C", direction = -1)

# Gene expression
# Krt8 (Fig.2c)
FeaturePlot(seuset, features="Krt8", split.by="sample")
p & scale_colour_gradientn(colours = rev(myPalette))

# Esr1 (Fig. 2d)
p <- FeaturePlot(seuset, features="Esr1", split.by="sample", min.cutoff = "q10", max.cutoff = "q90")
p & scale_colour_gradientn(colours = rev(myPalette))

# Kit (Fig. 2e)
p <- FeaturePlot(seuset, features="Kit", split.by="sample", min.cutoff = "q10", max.cutoff = "q95")
p & scale_colour_gradientn(colours = rev(myPalette))

# Itgb3 (Fig.2f)
p <- FeaturePlot(seuset, features="Itgb3", split.by="sample", min.cutoff = "q10", max.cutoff = "q60")
p & scale_colour_gradientn(colours = rev(myPalette))

# Krt14 (Fig.2g)
FeaturePlot(seuset, features="Krt14", split.by="sample")
p & scale_colour_gradientn(colours = rev(myPalette))

# Mki67 (Fig.2h)
FeaturePlot(seuset, features="Mki67", split.by="sample")
p & scale_colour_gradientn(colours = rev(myPalette))