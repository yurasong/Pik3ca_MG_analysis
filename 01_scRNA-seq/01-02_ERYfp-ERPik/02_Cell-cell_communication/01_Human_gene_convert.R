#####################################################################################
# Convert mouse genes to human orthologs and rebuild Seurat object
#
# Description:
#   Converts murine gene symbols to human orthologs using homologene mapping,
#   rebuilds a Seurat object with human gene names, and performs standard
#   pre-processing including normalization, feature selection, scaling,
#   dimensional reduction, clustering, and UMAP visualization.
#   The conversion is performed prior to downstream analyses to enable
#   cross-species comparison or integration with human datasets.
#
# Workflow overview:
#     1. Load annotated mouse Seurat object
#     2. Map mouse genes to human orthologs using homologene
#     3. Subset count matrix to mapped genes and rename features
#     4. Reconstruct Seurat object with identical metadata
#     5. Run standard Seurat preprocessing and clustering pipeline
#
# Inputs:
#   - ERPik_Stroma_annotate.rds: Annotated murine stromal Seurat object (ER Pik Yfp dataset)
#
# Outputs:
#   - ERPik_Stroma_human_converted.rds: Human-ortholog-converted Seurat object with UMAP and clusters
#
#
# Dependencies:
#   Seurat, tidyverse, dplyr, ggplot2, Matrix,
#   biovizBase, biomaRt, homologene
#
# Notes:
#   - Only one-to-one ortholog mappings are retained
#   - Duplicate mouse or human gene mappings are removed
#   - Metadata from original object is preserved exactly
#####################################################################################

# Library

library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(biovizBase)
library(Matrix)
library(biomaRt)
library(homologene)

# Data

erpik <- readRDS("ERPik_Stroma_annotate.rds")

# Convert to human genes

mouse_genes <- rownames(GetAssayData(erpik, assay = "RNA", layer = "data"))

map_df <- mouse2human(mouse_genes) # Get the homologene map locally
colnames(map_df) <- c("mouse","human")
map_df <- map_df[!duplicated(map_df$mouse) & !duplicated(map_df$human), ]

cts <- GetAssayData(erpik, assay = "RNA", layer = "counts")
if (is.null(cts)) {
  cts <- GetAssayData(erpik, assay = "RNA", layer = "data")
}

## Subset to only the mouse genes you mapped
keep   <- intersect(rownames(cts), map_df$mouse)
cts_sub <- cts[keep, , drop = FALSE]

## Rename the rows all at once
new_names <- setNames(map_df$human, map_df$mouse)
rownames(cts_sub) <- new_names[ rownames(cts_sub) ]

# Build a fresh Seurat object with the exact same metadata
seurat_human <- CreateSeuratObject(
  counts    = cts_sub,
  meta.data = erpik@meta.data,
  project   = "Mouse_to_Human"
)

## Proceed human-converted data

seurat_human <- NormalizeData(seurat_human, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_human <- FindVariableFeatures(seurat_human, selection.method = "vst", nfeatures = 2000)
seurat_human <- ScaleData(seurat_human)
seurat_human <- RunPCA(seurat_human, features = VariableFeatures(object = seurat_human))

seurat_human <- FindNeighbors(seurat_human, dims = 1:40)
seurat_human_res0p5 <- FindClusters(seurat_human, resolution = 0.5)
seurat_human_res0p5 <- RunUMAP(seurat_human_res0p5, dims = 1:40)

seurat_human_res0p5@active.ident =factor(as.character(seurat_human_res0p5@meta.data$cell_type))
names(seurat_human_res0p5@active.ident) = rownames(seurat_human_res0p5@meta.data)

# Export object

saveRDS(seurat_human_res0p5, "ERPik_Stroma_human_converted.rds")