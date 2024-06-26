# Since the data has been dealt with Signac, we have to create anndata (.h5ad) format for transcriptome.
# You need expression matrix and metadata including cellID and cell_type (or at least seurat clusters).

import pandas as pd
import scanpy as sc

data = pd.read_csv("exprMat.tsv", sep='\t', header=0, index_col=0)
meta = pd.read_csv("metadata_multiome.tsv", sep='\t', header=0, index_col=0)

adata = sc.AnnData(X=data, obs=meta)
adata.write("RNA.h5ad")
