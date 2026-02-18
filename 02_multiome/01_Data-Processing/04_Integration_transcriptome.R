#####################################################################################
# File: 04_integrate_transcriptome.R
#
# Description:
#   1. Load three pre-filtered multiome Seurat objects (K8Pik, K8Pik_Klf5KO, CTL)
#   2. Merge and integrate RNA assays across samples using Seurat integration workflow
#   3. Perform PCA, UMAP, neighbor graph, and clustering at resolution 0.5
#   4. Export cluster UMAPs, per-sample UMAPs, metadata, and top marker genes
#   5. Score cell cycle phases and plot Phase UMAP
#   6. Plot canonical marker gene expression (violin plots)
#   7. Sub-cluster specified populations and recluster at resolution 0.3
#   8. Merge subcluster metadata back into main object, update cluster labels
#   9. Rename idents to cell-type labels, generate final UMAPs (with and without legend, split by sample)
#  10. Save final annotated Seurat object and metadata CSV
#
# Inputs: Individual files which are proceeded
#   - K8pik_multiome.rds
#   - K8pik_klf5ko_multiome.rds
#   - ctl_multiome.rds
#
# Dependencies:
#   Signac, Seurat, tidyverse, dplyr, ggplot2,
#   EnsDb.Mmusculus.v79, BSgenome.Mmusculus.UCSC.mm10,
#   biovizBase, Matrix, DoubletFinder, forcats
#
#####################################################################################


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

# Dataset: Doublet filtered, no fibroblast contaminated (all QC done upfront)
k8pik <- readRDS("K8pik_multiome.rds")
klf5ko <- readRDS("K8pik_klf5ko_multiome.rds")
ctl <- readRDS("ctl_multiome.rds")

DefaultAssay(k8pik) <- "RNA"
DefaultAssay(klf5ko) <- "RNA"
DefaultAssay(ctl) <- "RNA"

# Set-up for integration
k8pik$sample <- "K8Pik"
klf5ko$sample <- "K8Pik_Klf5KO"
ctl$sample <- "CTL"

tmp <- merge(k8pik, klf5ko, add.cell.ids=c("K8Pik", "K8Pik_Klf5KO"))
fullset <- merge(tmp, ctl, project = "Multiome")

## Remove all intermediate dataset
rm(tmp)
rm(k8pik)
rm(klf5ko)
rm(ctl)

# Integration with Seurat

## Split the datasets into a list of Seurat objects:
multiome.list <- SplitObject(fullset, split.by = "sample")

multiome.list <- lapply(X = multiome.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

## select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = multiome.list)

## Integration
multiome.anchors <- FindIntegrationAnchors(object.list = multiome.list, anchor.features = features) 
multiome.combined <- IntegrateData(anchorset = multiome.anchors)

## Integrated analysis

DefaultAssay(multiome.combined) <- "integrated"

multiome.combined <- ScaleData(multiome.combined, verbose = FALSE)
multiome.combined <- RunPCA(multiome.combined, npcs = 30, verbose = FALSE)
ElbowPlot(multiome.combined)

multiome.combined <- RunUMAP(multiome.combined, reduction = "pca", dims = 1:30)
multiome.combined <- FindNeighbors(multiome.combined, reduction = "pca", dims = 1:30)

# Application of resolution

dir.create("res0p5") # To export plots

multiome.combined5 <- FindClusters(multiome.combined, resolution = 0.5)

DimPlot(multiome.combined5, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("res0p4/UMAP_Integrated_multiome_RNA_res0p5.pdf", height=9, width=12, units="cm")

DimPlot(multiome.combined5, reduction = "umap", split.by = "sample")
ggsave("res0p4/UMAP_Integrated_res0p5_by_sample.pdf", height=9, width=27, units="cm")

## Export meta.data
write.table(multiome.combined5@meta.data, "res0p4/metadata_Integrated_multiome_RNA_res0p5.csv", quote=F, sep="\t", row.names=F, col.names=T)

## Marker genes per cluster
markers <- FindAllMarkers(multiome.combined5, logfc.threshold = 0.5, min.pct = 0.25, only.pos = T)

markers <- markers %>% dplyr::filter(p_val_adj < 0.01)
markers <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.table(markers, "Marker_Top100_integrated_RNA_res0p5.csv", quote=F, sep="\t", row.names=F, col.names=T)

# Cell cycle analysis
s.genes <-str_to_title(cc.genes$s.genes)
g2m.genes <-str_to_title(cc.genes$g2m.genes)

DefaultAssay(multiome.combined5)  <- "RNA"

cc_inte <- CellCycleScoring(multiome.combined5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
cc_inte <- ScaleData(cc_inte)
cc_inte <- RunPCA(cc_inte, features = c(s.genes, g2m.genes))

p <- DimPlot(cc_inte, label=T, group.by = "Phase", reduction = "umap")
ggsave("Cell_cycle.pdf", p, width=12, height=9, units="cm")

p

# Canonical Makrer Gene Expression
bc <- c("Krt5", "Krt14", "Trp63", "Acta2", "Myh11", "Mylk")
lcp <- c("Foxa1", "Prlr", "Pgr", "Esr1", "Prom1")
lcm <- c("Elf5", "Kit", "Cd14", "Itga2", "Csn3")

DefaultAssay(multiome.combined5) <- "RNA"

VlnPlot(multiome.combined5, features=bc, ncol=3, pt.size=0)
VlnPlot(multiome.combined5, features=lcp, ncol=3, pt.size=0)
VlnPlot(multiome.combined5, features=lcm, ncol=3, pt.size=0)

VlnPlot(multiome.combined5, features="Krt8", pt.size=0) # Pan-LC marker
VlnPlot(multiome.combined5, features=c("Mki67", "Top2a", "Cenpa"), pt.size=0) # Proliferation

# Sub-clustering
## Here we did sub-clustering of certain clusters which was mixture of hybrid and distinct cells.

DefaultAssay(multiome.combined5) <- "integrated"
subc <- subset(multiome.combined5, idents=c(4, 7, 10))

# Subclustered datasets proceeding

subc <- ScaleData(subc, verbose = FALSE)
subc <- RunPCA(subc, npcs = 30, verbose = FALSE)
subc <- RunUMAP(subc, reduction = "pca", dims = 1:30)

subc <- FindNeighbors(subc, reduction = "pca", dims = 1:30)

## Application of resolution

dir.create("Subc_res0p3")
subc3 <- FindClusters(subc, resolution = 0.3)

DimPlot(subc3, reduction = "umap", label = TRUE, repel = TRUE)
ggsave("Subc_res0p3/UMAP_Integrated_multiome_RNA_subc_res0p3.pdf", height=9, width=12, units="cm")

DimPlot(subc3, reduction = "umap", split.by = "sample")
ggsave("Subc_res0p3/UMAP_Integrated_multiome_RNA_subc_res0p3.pdf", height=9, width=27, units="cm")

## Gene expression check
DefaultAssay(subc3) <- "RNA"

FeaturePlot(subc3, features=c(bc, lcp, lcm, "Mki67", "Top2a", "Cenpa"), ncol = 5, min.cutoff = "q5")
VlnPlot(subc3, features=bc, ncol=3, pt.size=0)
VlnPlot(subc3, features=lcp, ncol=3, pt.size=0)
VlnPlot(subc3, features=lcm, ncol=3, pt.size=0)
VlnPlot(subc3, features=c("Mki67", "Top2a", "Cenpa"), pt.size=0)

# Projection to original object

all_meta <- multiome.combined5@meta.data
subc_meta <- subc3@meta.data

subc_meta <- subc_meta[, c(1, 23)]
names(subc_meta) <- c("orig.ident", "subc_cluster")
names(subc_meta)

merged_meta <- merge(all_meta, subc_meta, by=0, all.x=T)
rownames(merged_meta) <- merged_meta$Row.names
merged_meta <- merged_meta[, -1]
names(merged_meta)

merged_meta <- merged_meta[, c(1:13, 18, 22, 24)]
names(merged_meta)

meta <- merged_meta %>%
  mutate(
    cluster_updated = if_else(
      !is.na(subc_cluster) & subc_cluster != "",
      paste0("s", subc_cluster),
      as.character(integrated_snn_res.0.5)
    )
  )

multiome.combined5@meta.data <- meta
names(meta)

multiome.combined5@active.ident =factor(as.character(multiome.combined5@meta.data$cluster_updated))
names(multiome.combined5@active.ident) = rownames(multiome.combined5@meta.data)

DefaultAssay(multiome.combined5) <- "integrated"

# Annotation

multiome.combined5 <- RenameIdents(multiome.combined5, 
                       `0` = "LC_ER-", 
                       `1` = "LC_ER+", 
                       `2` = "BCs",
                       `3` = "LC_ER+", 
                       `5` = "LC_ER-", 
                       `6` = "LC_ER+", 
                       `8` = "LC_ER+",
                     `9` = "Prolif",
                     `11` = "LC_ER+",
                     `12` = "BCs",
                     `s0` = "HY_ER+/ER-",
                       `s1` = "Myoepith", 
                       `s2` = "HY_BC/ER-",
                       `s3` = "HY_BC/ER+", 
                     `s4` = "HY_BC/ER+", 
                       `s5` = "HY_BC/ER+")

multiome.combined5$cell_type <- fct_infreq(multiome.combined5@active.ident)

multiome.combined5@active.ident =factor(as.character(multiome.combined5@meta.data$cell_type))
names(multiome.combined5@active.ident) = rownames(multiome.combined5@meta.data)

DimPlot(multiome.combined5, reduction = "umap", label = TRUE, repel = TRUE,
        cols=c("LC_ER-" = "#00BFC4",
               "BCs" = "#00BE67",
               "LC_ER+" = "#00A9FF",
               "HY_BC/ER-" = "#F8766D",
               "HY_ER+/ER-" = "#7CAE00",
               "Prolif" = "#FF61CC", 
               "Myoepith" = "#C77CFF",
               "HY_BC/ER+" = "#CD9600")) + NoLegend()
ggsave("UMAP.pdf", width = 12, height=9, units="cm")

DimPlot(multiome.combined5, label=F, 
        cols=c("LC_ER-" = "#00BFC4",
               "BCs" = "#00BE67",
               "LC_ER+" = "#00A9FF",
               "HY_BC/ER-" = "#F8766D",
               "HY_ER+/ER-" = "#7CAE00",
               "Prolif" = "#FF61CC", 
               "Myoepith" = "#C77CFF",
               "HY_BC/ER+" = "#CD9600")) + NoLegend() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            legend.position = "none",
            axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave("UMAP_noLegend.pdf", width = 10, height = 10, units = "cm")

DimPlot(multiome.combined5, reduction = "umap", split.by="sample",
        cols=c("LC_ER-" = "#00BFC4",
               "BCs" = "#00BE67",
               "LC_ER+" = "#00A9FF",
               "HY_BC/ER-" = "#F8766D",
               "HY_ER+/ER-" = "#7CAE00",
               "Prolif" = "#FF61CC", 
               "Myoepith" = "#C77CFF",
               "HY_BC/ER+" = "#CD9600")) + NoLegend()
ggsave("UMAP_split_by_sample.pdf", width = 27, height=9, units="cm")

saveRDS(multiome.combined5, "Multiome_integrated_annot.rds")
write.table(multiome.combined5@meta.data, "metadata.csv", quote=F, sep=",", row.names=T, col.names=T)