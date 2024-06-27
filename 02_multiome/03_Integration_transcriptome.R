# Library
library(Signac)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(biovizBase)
library(Matrix)

# Data

k8pik <- readRDS(".Annot_K8Pik_noFibro.rds")
klf5ko <- readRDS("Annot_K8Pik_Klf5KO_noFibro.rds")

DefaultAssay(k8pik) <- "RNA"
DefaultAssay(klf5ko) <- "RNA"

## Set-up for integration
k8pik$sample <- "K8Pik"
klf5ko$sample <- "K8Pik_Klf5KO"

fullset <- merge(k8pik, klf5ko, add.cell.ids=c("K8Pik", "K8Pik_Klf5KO"), project = "Multiome")

# Integration with Seurat
multiome.list <- SplitObject(fullset, split.by = "sample")

multiome.list <- lapply(X = multiome.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = multiome.list)

## Integration

multiome.anchors <- FindIntegrationAnchors(object.list = multiome.list, anchor.features = features)
multiome.combined <- IntegrateData(anchorset = multiome.anchors)

## Clustering

DefaultAssay(multiome.combined) <- "integrated"

multiome.combined <- ScaleData(multiome.combined, verbose = FALSE)
multiome.combined <- RunPCA(multiome.combined, npcs = 30, verbose = FALSE)
ElbowPlot(multiome.combined)

multiome.combined <- RunUMAP(multiome.combined, reduction = "pca", dims = 1:30)
multiome.combined <- FindNeighbors(multiome.combined, reduction = "pca", dims = 1:30)

# Resolution application

dir.create("res0p5")

multiome.combined5 <- FindClusters(multiome.combined, resolution = 0.5)
DimPlot(multiome.combined5, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(multiome.combined5, reduction = "umap", split.by = "sample")

# Cell cycle analysis

s.genes <-str_to_title(cc.genes$s.genes)
g2m.genes <-str_to_title(cc.genes$g2m.genes)

DefaultAssay(multiome.combined5)  <- "RNA"

cc_inte <- CellCycleScoring(multiome.combined5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cc_inte <- ScaleData(cc_inte)
cc_inte <- RunPCA(cc_inte, features = c(s.genes, g2m.genes))

p <- DimPlot(cc_inte, label=T, group.by = "Phase", reduction = "umap")
p

# Annotation

## Gene Expression
bc <- c("Krt5", "Krt14", "Trp63", "Acta2", "Myh11", "Mylk")
lcp <- c("Foxa1", "Prlr", "Pgr", "Esr1", "Prom1")
lcm <- c("Elf5", "Kit", "Cd14", "Itga2", "Csn3")
oth <- c("Krt8", "Mki67", "Klf5", "Krt6a", "Il33")

DefaultAssay(multiome.combined5) <- "RNA"

FeaturePlot(multiome.combined5, features=bc, ncol=3)
FeaturePlot(multiome.combined5, features=lcp, ncol=3)
FeaturePlot(multiome.combined5, features=lcm, ncol=3)
FeaturePlot(multiome.combined5, features=oth, ncol=3)

## Annotation
seuset <- RenameIdents(multiome.combined5, 
                       `0` = "LC_ER-", 
                       `1` = "Immature_BC", 
                       `2` = "LC_ER+",
                       `3` = "HY_ER+/ER-", 
                       `4` = "HY_ER+/BC", 
                       `5` = "Prolif", 
                       `6` = "Myoepith", 
                       `7` = "HY_ER+/BC", 
                       `8` = "LC_ER+_Klf5_specific",
                       `9` = "HY_ER+/BC", 
                       `10` = "LC_ER-", 
                       `11` = "HY BC/ER-",
                       `12` = "LC_ER+")

seuset$cell_type <- fct_infreq(seuset@active.ident)

DimPlot(seuset, label=T,
        cols=c("LC_ER-" = "#00BFC4",
               "Immature_BC" = "#00BE67",
               "LC_ER+" = "#00A9FF",
               "HY BC/ER-" = "#F8766D",
               "HY_ER+/ER-" = "#7CAE00",
               "Prolif" = "#FF61CC",
               "Myoepith" = "#C77CFF",
               "HY_ER+/BC" = "#CD9600",
               "LC_ER+_Klf5_specific" = "#800080")) + NoLegend() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("UMAP_RNA_Integrated_Label.pdf", width = 10, height = 10, units = "cm")

DimPlot(seuset, label=F, split.by="sample",
        cols=c("LC_ER-" = "#00BFC4",
               "Immature_BC" = "#00BE67",
               "LC_ER+" = "#00A9FF",
               "HY BC/ER-" = "#F8766D",
               "HY_ER+/ER-" = "#7CAE00",
               "Prolif" = "#FF61CC",
               "Myoepith" = "#C77CFF",
               "HY_ER+/BC" = "#CD9600",
               "LC_ER+_Klf5_specific" = "#800080")) + NoLegend() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("UMAP_Integrated-split.pdf", width = 20, height = 10, units = "cm")


# Save data
saveRDS(seuset, "Annot_multiome_integrated.rds")

