#####################################################################################
# File: 00_anndata_create.py
#
# Description:
#   Converts Seurat‐derived expression matrix and metadata into an AnnData
#   object and writes it to an H5AD file for downstream Scanpy analyses.
#
# Usage:
#   python 08_create_h5ad.py
#
# Inputs: Retrieve it from seurat object before running this workflow.
#   - exprMat.tsv               : Tab‐delimited gene expression matrix
#                                 (genes × cells, header row with cell IDs,
#                                 first column with gene names)
#   - metadata_multiome.tsv     : Tab‐delimited cell metadata
#                                 (rows keyed by cell IDs, includes cell_type or cluster)
#
# Outputs:
#   - RNA.h5ad                  : AnnData object saved in H5AD format
#
# Dependencies:
#   pandas, scanpy
#####################################################################################

import pandas as pd
import scanpy as sc

# Read expression matrix (genes × cells)
data = pd.read_csv("exprMat.tsv", sep="\t", header=0, index_col=0)

# Read cell metadata (cells × annotations)
meta = pd.read_csv("metadata_multiome.tsv", sep="\t", header=0, index_col=0)

# Create AnnData object
adata = sc.AnnData(X=data, obs=meta)

# Write to H5AD
adata.write("RNA.h5ad")