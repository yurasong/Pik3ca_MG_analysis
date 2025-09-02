#####################################################################################
# 05_03_01_input_file_from_Seurat_object.R
# Here, we used data object which slingshot pseudotime is embedded.
# Prepare input files for scanpy PAGA
#####################################################################################

# Library
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

seuset <- readRDS("Slingshot_embedded_integrated.rds")

# Gene expression matrix
exprMat <- as.matrix(seuset@assays$RNA@data)
write.table(data.frame("cell_id"=rownames(exprMat), exprMat), "ExprMat_row_gene.tsv", row.names=F, col.names=T, sep="\t", quote=F)

# Metadata
meta <- seuset@meta.data
write.table(meta, "metadata_including_slingshot_Pseudotime.csv", quote=F, sep=",", row.names=T, col.names=T)

# UMAP coorinate
umap_coord <- seuset@reductions$umap@cell.embeddings
write.table(umap_coord, "UMAP_coordinate.csv", quote=F, sep=",", row.names=T, col.names=T)