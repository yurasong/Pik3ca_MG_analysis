# System setting

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

import pickle
infile = open('/globalscratch/ysong/K8Pik/eGRN_inference/K8Pik_scplus_obj.pkl', 'rb')
scplus_obj = pickle.load(infile)
infile.close()

## Enhancer-to-gene model

from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, GBM_KWARGS

calculate_regions_to_genes_relationships(scplus_obj,
                    ray_n_cpu = 10,
                    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
                    importance_scoring_method = 'GBM',
                    importance_scoring_kwargs = GBM_KWARGS)

import pickle
pickle.dump(scplus_obj, open(os.path.join(work_dir, 'eGRN_inference/K8Pik_scplus_obj.pkl'), 'wb'))

# Infer TF to gene relationship

from scenicplus.TF_to_gene import *
tf_file = '/home/ulb/iribhm/ysong/database_scenicplus/allTFs_mm.txt'
calculate_TFs_to_genes_relationships(scplus_obj,
                    tf_file = tf_file,
                    ray_n_cpu = 10,
                    method = 'GBM',
                    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
                    key= 'TF2G_adj')

## Save
import pickle
pickle.dump(scplus_obj, open(os.path.join(work_dir, 'eGRN_inference/K8Pik_scplus_obj.pkl'), 'wb'))

# Build eGRNs

# Load functions
from scenicplus.grn_builder.gsea_approach import build_grn

build_grn(scplus_obj,
         min_target_genes = 10,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.05,
         keep_extended_motif_annot = True,
         merge_eRegulons = True,
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = 10,
         _temp_dir = os.path.join(tmp_dir, 'ray_spill'))

## To access the eGRNs:

import dill
with open('/globalscratch/ysong/K8Pik/eGRN_inference/K8Pik_scplus_obj.pkl', 'wb') as f:
  dill.dump(scplus_obj, f)



