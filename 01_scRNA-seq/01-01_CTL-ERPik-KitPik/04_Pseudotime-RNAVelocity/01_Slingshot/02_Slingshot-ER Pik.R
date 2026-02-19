#####################################################################################
# Slingshot trajectory inference on ER_Pik epithelial cells and pseudotime-GAM dynamics
#
# Description:
#   Runs Slingshot trajectory inference on ER_Pik epithelial cells (non-proliferative)
#   using PCA embeddings from an integrated Seurat object. The workflow converts
#   Seurat data into a SingleCellExperiment, fits Slingshot lineages under both
#   free and fixed topology settings, and stores inferred pseudotime back into
#   Seurat metadata. Subsequently, gene expression dynamics along pseudotime are
#   modeled using a GAM-based differential expression framework on variable genes
#   (custom functions from sourced scripts), and results are stored in the Seurat
#   object for downstream plotting and interpretation.
#
# Workflow overview:
#     1. Load integrated epithelial Seurat object
#     2. Subset to ER_Pik sample and remove proliferating cells ("Prol")
#     3. Create output directories for figures and results
#     4. Source custom pseudotime analysis and plotting scripts
#     5. Build SingleCellExperiment with:
#          - log-normalized expression (Seurat "data" slot)
#          - PCA coordinates (PC1–PC30)
#          - metadata (seurat_clusters, cell_type)
#     6. Quality control: visualize PCA colored by cell_type and marker genes
#     7. Run Slingshot:
#          a) free topology (unspecified start/end)
#          b) fixed topology with defined start/end clusters
#     8. Remove any previous slingshot pseudotime columns and add new ones to Seurat
#     9. Fit pseudotime-dependent gene expression using GAM along lineage 1
#        on Seurat variable features, and store DE results in `seuset@misc`
#
# Inputs:
#   - Epithelial_Integrated_Annot.rds
#       Integrated epithelial Seurat object containing PCA reduction, variable
#       features, and metadata columns including `cell_type` and `sample`
#   - scripts/pseudotemporalDynamics.R
#       Custom functions for pseudotime GAM fitting (e.g., fitPseudotimeGAM_nopal)
#   - scripts/plots.R
#       Custom plotting utilities used downstream
#
# Outputs:
#   - Slingshot pseudotime annotations added to Seurat meta.data:
#       slingPseudotime_1 (and additional lineages if present)
#   - PCA diagnostic plots (cell_type and marker overlays) for embedding sanity check
#   - Stored pseudotime-DE results:
#       seuset@misc$PT_DE_results$slingshot
#         • TI_method = "Slingshot"
#         • DE_method = "GAM"
#         • results$lineage1 = gam_fdr_t1
#
#
# Dependencies:
#   Seurat, SingleCellExperiment, scater, slingshot, tidyverse,
#   ggpubr, ggthemes, GGally, RColorBrewer, gam
#
# Notes:
#   - Slingshot is run on PCA space (not UMAP) for trajectory fitting stability.
#####################################################################################


# Preparation
## Library
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(GGally)
library(RColorBrewer)
library(gam, quietly = TRUE)

## Data
seuset <- readRDS("Epithelial_Integrated_Annot.rds")

## Filter out proliferating cells

seuset <- seuset[, seuset$cell_type != "Prol"]
seuset <- seuset[, seuset$sample == "ER_Pik"]

## Configuration
fig_path <- file.path("figures_ER_Pik")
dir.create(fig_path, recursive = TRUE)
res_path <- file.path("results_ER_Pik")
dir.create(res_path, recursive = TRUE)

sling_fig_path <- file.path(fig_path, "slingshot-PCA_ER_Pik")
dir.create(sling_fig_path)

## Load script
source("scripts/pseudotemporalDynamics.R")
source("scripts/plots.R")

# Slingshot on PCA

## Generate input data
exprs_data <- GetAssayData(seuset, slot = "data")
meta <- seuset@meta.data[c("seurat_clusters", "cell_type")]
pca <- Embeddings(seuset, reduction = "pca")[, 1:30]

sce <- SingleCellExperiment(
  assays = exprs_data,
  colData = meta
)

logcounts(sce) <- exprs_data
reducedDim(sce, "PCA") <- pca

### To check whether the cell types are correER_Piky embedded
p1 <- plotPCA(sce, colour_by = "cell_type")
p2 <- plotPCA(sce, colour_by = "Krt14")
p3 <- plotPCA(sce, colour_by = "Acta2")
p4 <- plotPCA(sce, colour_by = "Prlr")
CombinePlots(list(p1, p2, p3, p4), ncol = 2)

## Apply Slingshot with free topology
sling_free <- slingshot(sce, clusterLabels = "cell_type", reducedDim = "PCA")
SlingshotDataSet(sling_free)

## Apply Slingshot with Fixed topology
sling_fixed <- slingshot(sce, clusterLabels = "cell_type", reducedDim = "PCA", 
                         start.clus = c("ER+ LCs Foxa1", "ER+ LCs Sca1"), end.clus = "Myo")
SlingshotDataSet(sling_fixed)

# Embedding
n_lineages <- SlingshotDataSet(sling_fixed)@lineages %>% length()

## remove old slingshot results
seuset@meta.data[which(
  str_detect(names(seuset@meta.data), 
             "^slingPseudotime"))] <- NULL

## add new results
sling_names <- names(sling_fixed@colData) %>%
  .[str_detect(., "^slingPseudotime")]

seuset@meta.data[sling_names] <- as.data.frame(sling_fixed@colData[sling_names])

# Pseudotemporal gene dynamics
## set pseudotime variables
t1 <- seuset$slingPseudotime_1

## extract gene expression data for var.genes
Y <- FetchData(seuset, vars = VariableFeatures(seuset), slot = "data")
gam_fdr_t1 <- fitPseudotimeGAM_nopal(t1, Y) 

## Add on Seurat object
seuset@misc$PT_DE_results$slingshot <- list(
  TI_method = "Slingshot",
  DE_method = "GAM",
  results = list(lineage1 = gam_fdr_t1)
  )