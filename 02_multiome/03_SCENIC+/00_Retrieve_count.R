#####################################################################################
# File: 00_Retrieve_count.R
#
# Description:
#   Export cell metadata, gene expression matrix, and ATAC peak counts
#   from a Seurat object into TSV formats for pycisTopic input.
#
# Inputs:
#   - seurat_object.rds : Processed and annotated Seurat multiome object 
#                        (must contain RNA and ATAC assays)
#
# Outputs:
#   - metadata.tsv      : Cell metadata table (rows = cell IDs)
#   - exprMat.tsv       : Tab-delimited gene expression matrix 
#                        (first column = cell_id, remaining columns = gene counts)
#   - CellID.txt        : List of ATAC cell IDs with “:” as separator
#   - count_table.tsv   : Peak-by-cell count matrix with matching CellIDs
#
# Dependencies:
#   Signac, Seurat, dplyr, tidyverse, magrittr, data.table,
#   patchwork, clustree
#####################################################################################



library(Signac)
library(Seurat)
library(dplyr)
library(patchwork)
library(clustree)
library(tidyverse)
library(magrittr)
library(data.table)

seuset <- readRDS("seurat_object.rds") # Swap with object which needs to be used.

# Metadata
write.table(seuset@meta.data, "metadata.tsv", quote=F, sep="\t", row.names=T, col.names=T)

# Gene expression
exprMat <- t(as.matrix(seurat_object@assays$RNA@data))
write.table(data.frame("cell_id"=rownames(exprMat), exprMat), "exprMat.tsv", row.names=F, col.names=T, sep="\t", quote=F)

# ATAC-seq
cellIDs <- rownames(seuset@assays$peaks@counts)
cellIDs_fixed <- sub("-", ":", cellIDs)
write.table(cellIDs_fixed,
            file     = "CellID.txt",
            quote    = FALSE,
            col.names = FALSE,
            row.names = FALSE)

updated <- cbind(cellIDs_fixed, exprMat)
write.table(updated, "count_table.tsv", row.names=F, col.names=T, sep="\t", quote=F)