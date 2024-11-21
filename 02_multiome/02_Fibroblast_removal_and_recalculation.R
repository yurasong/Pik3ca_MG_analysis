#####################################################################################
# This code is applied across the data.
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

# Data

seuset <- readRDS("01_Multiome_K8Pik.rds") # Object which already passed QC

DefaultAssay(seuset) <- "RNA"
FeaturePlot(seuset, features=c("Vim", "Fn1", "Pdgfra"), ncol=3, pt.size=0.8)

## Subsetting
subc <- subset(seuset, idents=c(0, 1, 2, 3, 4, 5, 6, 7)) # Take the clusters which does not express fibroblast markers

# subclustered datasets proceeding

## Gene expression

subc <- FindVariableFeatures(subc, selection.method = "vst", nfeatures = 1000)

all.genes <- rownames(subc)
subc <- ScaleData(subc, features = all.genes)

subc <- RunPCA(subc, features = VariableFeatures(object = subc))

ElbowPlot(subc,ndim=50) # For diagnosis
subc <- FindNeighbors(subc, dims = 1:40, k.param = 5)

## DNA accessibility data processing

DefaultAssay(subc) <- "peaks"
subc <- FindTopFeatures(subc, min.cutoff = 5)
subc <- RunTFIDF(subc)
subc <- RunSVD(subc) 

DepthCor(subc)

# Joint UMAP calculation

## build a joint neighbor graph using both assays
subc <- FindMultiModalNeighbors(
  object = subc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:40), 
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

DefaultAssay(subc) <- "RNA"
subc <- FindClusters(subc, resolution = 0.6)

# UMAP visualization
subc <- RunUMAP(
  object = subc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(subc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

# Check annotation

## Gene expression

bc <- c("Krt5", "Krt14", "Trp63") # BC markers
myo <- c("Acta2", "Myh11", "Mylk") # Myoepithelial markers
lcp <- c("Foxa1", "Prlr", "Pgr", "Esr1", "Prom1") # ER+ LCs markers
lcm <- c("Elf5", "Kit", "Cd14", "Itga2", "Csn3") # ER- LCs markers
oth <- c("Krt8", "Mki67") # Pan LC, proliferating
fibro <- c("Vim", "Fn1", "Pdgfra", "Fbln2") # Fibroblast markers

DefaultAssay(seuset) <- "RNA"

FeaturePlot(subc, features=bc, ncol=3) 
FeaturePlot(subc, features=myo, ncol=3) 
FeaturePlot(subc, features=lcp, ncol=3)
FeaturePlot(subc, features=lcm, ncol=3)
FeaturePlot(subc, features=oth, ncol=2)
FeaturePlot(subc, features=fibro, ncol=2) # No clusters will express these markers.

# Annotation

subc <- RenameIdents(subc,
                     `0` = "Immature_BC",
                     `1` = "LC_ER+",
                     `2` = "LC_ER-",
                     `3` = "HY_ER/ER-",
                     `4` = "Myoepith",
                     `5` = "Prolif",
                     `6` = "LC_ER+",
                     `7` = "HY_BC/ER-",
                     `8` = "HY_BC/ER+",
                     `9` = "LC_ER+",
                     `10` = "HY_BC/ER-")
subc$cell_type <- fct_infreq(subc@active.ident)
subc@active.ident =factor(as.character(subc@meta.data$cell_type))
names(subc@active.ident) = rownames(subc@meta.data)

## UMAP with annotation

DimPlot(seuset, label=F,
        cols=c("LC_ER-" = "#00BFC4",
               "Immature_BC" = "#00BE67",
               "LC_ER+" = "#00A9FF",
               "HY_BC/ER-" = "#F8766D",
               "HY_ER/ER-" = "#7CAE00",
               "Prolif" = "#FF61CC",
               "Myoepith" = "#C77CFF",
               "HY_BC/ER+" = "#CD9600")) + 
  NoLegend() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave("UMAP_K8Pik.pdf", width = 10, height = 10, units = "cm")

# Marker genes for each cluster

markers <- FindAllMarkers(subc, logfc.threshold = 0.5, min.pct = 0.25, only.pos = T)
markers <- markers %>% dplyr::filter(p_val_adj < 0.01)
markers <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.table(markers, "Marker_Top100_significant_FC0p5_K8Pik_noFibro.csv", quote=F, sep="\t", row.names=F, col.names=T)

# Save object

saveRDS(subc, "Annot_K8Pik_noFibro.rds")

