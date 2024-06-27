#####################################################################################
# This code is applied on K8PIk.
# All parameters for QC and normalisation are exactly same for Klf5KO data.
# The reference of code : https://stuartlab.org/signac/articles/pbmc_multiomic.html
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
## load the RNA and ATAC data
counts <- Read10X_h5("filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5")
fragpath <- "ATACseq/atac_fragments.tsv.gz"

## Gene annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Assay definition
## create a Seurat object containing the RNA adata
seuset <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

## Add ATAC assay to the object
seuset[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

# Quality control

DefaultAssay(seuset) <- "ATAC"
seuset <- NucleosomeSignal(seuset)
seuset <- TSSEnrichment(seuset)

VlnPlot(
  object = seuset,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

## For detailed check

VlnPlot(seuset, features="nCount_RNA", pt.size = 0) + ylim(0, 50000) + NoLegend()
VlnPlot(seuset, features="nCount_RNA", pt.size = 0) + ylim(0, 5000) + NoLegend()
VlnPlot(seuset, features="nCount_ATAC", pt.size = 0) + ylim(0, 200000) + NoLegend()

## filter out low quality cells
seuset <- subset(
  x = seuset,
  subset = nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    nCount_RNA < 35000 &
    nCount_RNA > 2000 &
    nCount_ATAC < 100000 &
    nCount_ATAC > 1000
)

VlnPlot(
  object = seuset,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# Peak calling

## call peaks using MACS2
peaks <- CallPeaks(seuset)

## remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE) 

## quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(seuset),
  features = peaks,
  cells = colnames(seuset)
)

## create a new assay using the MACS2 peak set and add it to the Seurat object
seuset[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

# Data processing
## Gene expression data processing

DefaultAssay(seuset) <- "RNA"

seuset <- NormalizeData(seuset, normalization.method = "LogNormalize", scale.factor = 10000)
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(seuset)
seuset <- ScaleData(seuset, features = all.genes)

seuset <- RunPCA(seuset, features = VariableFeatures(object = seuset))

ElbowPlot(seuset,ndim=50) # for diagnosis and decide how many PCs to take

seuset <- FindNeighbors(seuset, dims=1:40)

## DNA accessibility data processing

DefaultAssay(seuset) <- "peaks"
seuset <- FindTopFeatures(seuset, min.cutoff = 5)
seuset <- RunTFIDF(seuset)
seuset <- RunSVD(seuset) 

DepthCor(seuset) # If the first LSI is 1 or -1, it indicates technical variant (strong correlation between the first LSI component and the total number of counts for the cell)

# Joint UMAP visualization

seuset <- FindMultiModalNeighbors(
  object = seuset,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:40, 2:40), 
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

seuset <- FindClusters(seuset, resolution = 0.5)

## Build UMAP

seuset <- RunUMAP(
  object = seuset,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(seuset, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

# Check annotation

## Gene expression

bc <- c("Krt5", "Krt14", "Trp63") # BC markers
myo <- c("Acta2", "Myh11", "Mylk") # Myoepithelial markers
lcp <- c("Foxa1", "Prlr", "Pgr", "Esr1", "Prom1") # ER+ LCs markers
lcm <- c("Elf5", "Kit", "Cd14", "Itga2", "Csn3") # ER- LCs markers
oth <- c("Krt8", "Mki67") # Pan LC, proliferating
fibro <- c("Vim", "Fn1", "Pdgfra", "Fbln2") # Fibroblast markers

DefaultAssay(seuset) <- "RNA"

FeaturePlot(seuset, features=bc, ncol=3) 
FeaturePlot(seuset, features=myo, ncol=3) 
FeaturePlot(seuset, features=lcp, ncol=3)
FeaturePlot(seuset, features=lcm, ncol=3)
FeaturePlot(seuset, features=oth, ncol=2)
FeaturePlot(seuset, features=fibro, ncol=2) # Need to remove the clusters which are expressing on this plots

## Linking peaks to genes

DefaultAssay(seuset) <- "peaks"

### compute the GC content for each peak
seuset <- RegionStats(seuset, genome = BSgenome.Mmusculus.UCSC.mm10)

seuset <- LinkPeaks(
  object = seuset,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Krt14", "Acta2",
                "Krt8", "Prlr", "Esr1", 
                "Kit", "Csn3")
)

CoveragePlot(
  object = seuset,
  region = "Krt14",
  features = "Krt14",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)

CoveragePlot(
  object = seuset,
  region = "Acta2",
  features = "Acta2",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)

CoveragePlot(
  object = seuset,
  region = "Krt8",
  features = "Krt8",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)

CoveragePlot(
  object = seuset,
  region = "Prlr",
  features = "Prlr",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)

CoveragePlot(
  object = seuset,
  region = "Esr1",
  features = "Esr1",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)

CoveragePlot(
  object = seuset,
  region = "Kit",
  features = "Kit",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)

CoveragePlot(
  object = seuset,
  region = "Csn3",
  features = "Csn3",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)


# Save object

saveRDS(seuset, "01_Multiome_K8Pik.rds")

