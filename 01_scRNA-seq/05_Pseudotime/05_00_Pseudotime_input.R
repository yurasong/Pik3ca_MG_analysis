#####################################################################################
# 05_00_Pseudotime_input.R
# This code is used for generating data object for pseudotime analysis.
#####################################################################################

# library
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)

seuset <- readRDS("Seurat_integrated_res0p5.rds")

seuset_sub <- subset(seuset, idents=c("Late_HY", "Early_HY", "Myo", "LC_ER+_Sca1", "LC_ER+_Foxa1"))

saveRDS(seuset_sub, "00_Pseudotime_input.rds")