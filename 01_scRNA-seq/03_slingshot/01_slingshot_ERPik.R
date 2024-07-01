#####################################################################################
# This code is applied on ERPik data.
#####################################################################################

# Library
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(destiny) # Manual installation required if your bioconductor version > 3.10
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(GGally)
library(RColorBrewer)
library(gam, quietly = TRUE)

theme_set(theme_light())

# Data
seuset <- readRDS("ERPik_annotated.rds")
seuset <- subset(seuset, idents=c("LC_ER+_Sca1", "LC_ER+_Foxa1", "BCs", 
                                  "Myoepith", "HY_BC/ER+", "LC_ER-"))

## Load the scripts
source("scripts/pseudotemporalDynamics.R")
source("scripts/plots.R")

# Slingshot on PCA

## Input data
exprs_data <- GetAssayData(seuset, slot = "data")
meta <- seuset@meta.data[c("seurat_clusters", "cell_type")]
pca <- Embeddings(seuset, reduction.type = "pca")[, 1:30]

sce <- SingleCellExperiment(
  assays = exprs_data,
  colData = meta)

logcounts(sce) <- exprs_data
reducedDim(sce, "PCA") <- pca

# Apply Slingshot

sling_fixed <- slingshot(sce, clusterLabels = "cell_type", reducedDim = "PCA", 
                         start.clus = c("LC_ER+_Sca1", "LC_ER+_Foxa1"), end.clus="Myoepith")

SlingshotDataSet(sling_fixed)

## Add slingshot pseudo-timelines to Seurat object and plot ordering of cells
n_lineages <- SlingshotDataSet(sling_fixed)@lineages %>% length()

### remove old slingshot results
seuset@meta.data[which(
  str_detect(names(seuset@meta.data), 
             "^slingPseudotime"))] <- NULL

### add new results
sling_names <- names(sling_fixed@colData) %>%
  .[str_detect(., "^slingPseudotime")]

seuset@meta.data[sling_names] <- as.data.frame(sling_fixed@colData[sling_names])

# Pseudotemporal gene dynamics
## set pseudotime variables
t1 <- seuset$slingPseudotime_1

## extract gene expression data for var.genes
Y <- FetchData(seuset, vars = VariableFeatures(seuset), slot = "data")

gam_fdr_t1 <- fitPseudotimeGAM_nopal(t1, Y) 

## add results to seurat
seuset@misc$PT_DE_results$slingshot <- list(
  TI_method = "Slingshot",
  DE_method = "GAM",
  results = list(lineage1 = gam_fdr_t1)
)

## Plotting on UMAP
sling_names <- names(seuset@meta.data) %>%
  .[str_detect(., "^slingPseudotime")]

sling_umaps <- FeaturePlot(seuset,
                           reduction = "umap",
                           features = sling_names, combine = FALSE
) %>%
  map(~ . + scale_color_viridis_c(option = "inferno", na.value = "light grey") +
        theme_void() + theme(aspect.ratio = 1))

CombinePlots(plots = sling_umaps, ncol = 1) %>% 
  annotate_figure(
    top = text_grob("Slingshot with fixed topology",
                    face = "bold", size = 16)
  )

## Plotting on PCA
sling_umaps <- FeaturePlot(seuset,
                           reduction = "pca",
                           features = sling_names, combine = FALSE
) %>%
  map(~ . + scale_color_viridis_c(option = "inferno", na.value = "light grey") +
        theme_void() + theme(aspect.ratio = 1))

CombinePlots(plots = sling_umaps, ncol = 1) %>% 
  annotate_figure(
    top = text_grob("Slingshot with fixed topology",
                    face = "bold", size = 16)
  )
