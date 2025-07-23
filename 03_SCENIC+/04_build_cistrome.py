#####################################################################################
# File: 04_build_cistrome.py
#
# Description:
#   1. Load scRNA-seq AnnData, cisTopic object, and motif enrichment results
#   2. Create a ScenicPlus object combining transcriptome and chromatin topic data
#   3. Generate and merge cistromes (eGRNs) for downstream regulatory network analysis
#   4. Save the ScenicPlus object for later eGRN inference steps
#
# Inputs:
#   - scRNA AnnData:      scRNA/RNA.h5ad
#   - cisTopic object:    output/K8Pik_cisTopicObject.pkl
#   - motif enrichment:   motifs/menr.pkl
#
# Outputs:
#   - ScenicPlus object:  eGRN_inference/K8Pik_scplus_obj.pkl
#   - Console log of merge_cistromes runtime
#
# Dependencies:
#   dill, scanpy, pandas, pyranges, scenicplus,
#   scenicplus.scenicplus_class.create_SCENICPLUS_object,
#   scenicplus.cistromes.merge_cistromes
#
# Notes:
#   - Ensure barcode mapping function matches your dataset conventions.
#   - Suppresses warnings for cleaner runtime logs.
#####################################################################################


import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
#_stderr = sys.stderr
#null = open(os.devnull,'wb')

work_dir = '/globalscratch/ysong/K8Pik'
tmp_dir = '/scratch/ysong'

# Data load

adata = sc.read_h5ad(os.path.join(work_dir, 'scRNA/RNA.h5ad'))
cistopic_obj = dill.load(open(os.path.join(work_dir, 'output/K8Pik_cisTopicObject.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np

# from scenicplus.scenicplus_class import create_SCENICPLUS_object
# help(create_SCENICPLUS_object)

scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata,
    cisTopic_obj = cistopic_obj,
    menr = menr,
    bc_transform_func = lambda x: f'{x}___cisTopic' 
#function to convert scATAC-seq barcodes to scRNA-seq ones
)

scplus_obj

# Generate cistrome

## Merge cistrome object

from scenicplus.cistromes import *
import time
start_time = time.time()
merge_cistromes(scplus_obj)
time = time.time()-start_time
print(time/60)

import pickle
pickle.dump(scplus_obj, open(os.path.join(work_dir, 'eGRN_inference/K8Pik_scplus_obj.pkl'), 'wb'))
