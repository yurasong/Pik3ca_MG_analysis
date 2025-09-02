#####################################################################################
# File: 02_AUC_LateHY.R
#
# Description:
#   1. Load processed Seurat object from Wu et al. (Nat Genet, 2021)
#   2. Build AUCell rankings and calculate AUC for the Late_HY_signature gene set
#   3. Add AUC assay to Seurat object and regress out feature count effects
#   4. Generate boxplots of adjusted AUC by patient and subtype
#
# Inputs:
#   - Wu_et_al_Cancer_epith_only.rds : Seurat object with RNA counts and metadata
#   - 01_Signature/Late_HY.txt      : List of genes defining the Late HY signature
#
# Outputs:
#   - Seurat object augmented with AUC assay (in memory)
#   - Boxplot of adjusted AUC by patient/subtype (plotted to current device)
#
# Dependencies:
#   dplyr, data.table, tidyverse, patchwork, Seurat,
#   AUCell, GSEABase, RColorBrewer, ggplot2, MASS
#####################################################################################

# Library
library(dplyr)
library(data.table)
library(tidyverse)
library(patchwork)
library(Seurat)

library(AUCell)
library(GSEABase)

library(RColorBrewer)
library(ggplot2)

library(MASS)

# Object
seuset <- readRDS("Wu_et_al_Cancer_epith_only.rds")

# Geneset allocation
early_hy <- read.delim("01_Signature/Late_HY.txt", h=F)
early_hy <- early_hy$V1
early_hy <- unique(early_hy)

geneSets <- GeneSet(early_hy, setName="Late_HY_signature")

# Running AUCell
exprMatrix <- as.matrix(seuset[["RNA"]]$counts)

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
response_AUC <- AUCell_calcAUC(geneSets, cells_rankings)

# Embedding

AUC_MA <- response_AUC@assays@data$AUC
AUC_MA <- as.matrix(AUC_MA)

aucs <- CreateAssayObject(counts = AUC_MA)

seuset[["AUC"]] <- aucs

# Correction by count
exprMatrix <- as.data.frame(seuset[["AUC"]]$counts)

meta <- seuset@meta.data[, c(1, 4, 8)]

merged_table <- merge(exprMatrix, meta, by=0)
rownames(merged_table) <- merged_table$Row.names
merged_table <- merged_table[, -1]

names(merged_table) <- c("AUC", "Patient", "nFeature_RNA", "subtype")

regressed <- rlm(AUC ~ nFeature_RNA, data=merged_table)

fitted <- mean(merged_table$AUC) + regressed$residuals[rownames(merged_table)]
merged_table[rownames(merged_table), "adjusted_AUC"] <- fitted

ggplot(merged_table, aes(x=Patient, y=adjusted_AUC, fill=subtype)) + 
  geom_boxplot() +
  scale_x_discrete(limits=c("CID3963", "CID4465", "CID4495", "CID44971", "CID44991", "CID4513",
                            "CID4515", "CID4523","CID3921", "CID4066", "CID45171", "CID3941", 
                            "CID3948", "CID4067", "CID4290A", "CID4461",
                            "CID4463", "CID4471", "CID4530N", "CID4535")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
