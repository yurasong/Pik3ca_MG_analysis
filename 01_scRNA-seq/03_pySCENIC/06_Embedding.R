#######################################################################################
# 06_embed_data_on_seurat.R
# This script will embed average AUC score into the Seurat object as a new assay.
#######################################################################################

# Preparations

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(dplyr)
library(patchwork)
library(clustree)
library(tidyverse)
library(pheatmap)

# Prepare the datasets

## Average result

scenic_results <- read.table("AUC_Average.csv", header=T, sep=",", row.names=1) # Avarage AUC value from 04_merge_results.py
colnames(scenic_results) <- sub("...$","",colnames(scenic_results))

## Seurat proceeded object

seuset <- readRDS("seurat_object.rds")

# Embedding on Seurat

# AUC matrix as the assay of seurat object
## To put the AUC matrix as additional assay on seurat object

seuset_auc <- seuset # To avoid issues on over-writing

auc_mtx <- t(scenic_results)
aucs <- CreateAssayObject(counts = auc_mtx)
seuset_auc[["AUC"]] <- aucs

DefaultAssay(seuset_auc) <- "AUC"