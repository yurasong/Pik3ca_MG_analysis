##########################################################################################################
# 05_01_Slinghost.R
# This code is used for inferring pseudotime trajectory on scRNA-seq data with slingshot R package.
##########################################################################################################

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
seuset <- readRDS("../00_Pseudotime_input.rds")

## Configuration
fig_path <- file.path("figures")
dir.create(fig_path, recursive = TRUE)
res_path <- file.path("results")
dir.create(res_path, recursive = TRUE)

sling_fig_path <- file.path(fig_path, "slingshot-PCA")
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

### To check whether the cell types are correctly embedded
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
                         start.clus = "LC_ER+_Foxa1", end.clus = "Myo")
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

p1 <- seuset@meta.data %>%
  dplyr::mutate(cell_type = fct_reorder(cell_type, slingPseudotime_1)) %>%
  ggplot(aes(slingPseudotime_1, cell_type, col = cell_type)) +
  geom_jitter() + 
  ggtitle("Slingshot cell ordering", subtitle = "Lineage 1") +
  scale_color_manual(values=c("#F08080", "#800000", "#BDB76B", "#32CD32", "#1E90FF"))

p1 # Will return ordering of cells by pseudotime calculated by Slingshot

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

# Export related files and embedded object
write.table(seuset@meta.data, "metadata_including_Pseudotime.csv", quote=F, sep=",", row.names=T, col.names = T)
saveRDS(seuset, "Slingshot_embedded_integrated.rds")