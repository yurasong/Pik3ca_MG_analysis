##########################################################################################
# File: 05_GetSearchSpace.sh
# Author: Yura Song
# Date: 2025-07-23
#
# Description:
#   Loads the ScenicPlus object and defines the enhancer-to-gene search space
#   using Ensembl Biomart (100 bp–100 kb upstream/downstream windows).
#   Must be run on the login node (internet access required).
#
# This script will execute:
#   - Load scplus_obj from eGRN_inference/K8Pik_scplus_obj.pkl
#   - Call get_search_space() with biomart_host, species, assembly, upstream/downstream
#   - Resave the updated scplus_obj to the same file
#
# Inputs:
#   - eGRN_inference/K8Pik_scplus_obj.pkl
#
# Outputs:
#   - Updated eGRN_inference/K8Pik_scplus_obj.pkl with search space defined
#
# Dependencies:
#   Python modules: dill, scanpy, numpy, pandas, pyranges, scenicplus.enhancer_to_gene
#
# Notes:
#   - Requires internet connection to Ensembl archive on login node.
##########################################################################################


# System setting

import dill
import scanpy as sc
import numpy as np
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

biomart_host = "http://sep2019.archive.ensembl.org/" # Most overlapping one (29476/32285)

# Since wrapper function does not work, we will return to differnt ways and use generating cistromes (https://scenicplus.readthedocs.io/en/latest/Scenicplus_step_by_step-RTD.html#3.-Infer-enhancer-to-gene-relationships).


scplus_obj = dill.load(open(os.path.join(work_dir, 'eGRN_inference/K8Pik_scplus_obj.pkl'), 'rb'))

from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS

# Get search space: Define search space, using 100kb up/down-stream the gene. 
# This step should be done on log-in node due to the internet connection
get_search_space(scplus_obj,
                 biomart_host = 'http://sep2019.archive.ensembl.org/',
                 species = 'mmusculus',
                 assembly = 'mm10',
                 upstream = [1000, 100000],
                 downstream = [1000, 100000])


import pickle
pickle.dump(scplus_obj, open(os.path.join(work_dir, 'eGRN_inference/K8Pik_scplus_obj.pkl'), 'wb'))
