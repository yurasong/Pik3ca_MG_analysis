seuset <- readRDS("seurat_object.RDS")

exprMat <- t(as.matrix(seuset@assays$RNA@data))
write.table(data.frame("cell_id"=rownames(exprMat), exprMat), "exprMat.tsv", row.names=F, col.names=T, sep="\t", quote=F)
