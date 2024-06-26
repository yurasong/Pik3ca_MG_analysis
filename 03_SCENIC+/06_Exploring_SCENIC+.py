# Exploring SCENIC+ results

import os
import warnings
warnings.filterwarnings("ignore")

# 1. Generate eRegulon metadata

import dill
infile = open('eGRN_inference/K8Pik_scplus_obj.pkl', 'rb')
scplus_obj = dill.load(infile)
infile.close()

from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')

## Take a look:
scplus_obj.uns['eRegulon_metadata'][0:10]

## Export metadata as table
import pandas as pd

meta_data = scplus_obj.uns['eRegulon_metadata']
meta_data.to_csv("eGRN_inference/K8Pik_eRegulon_metadata.csv", index=False)

# 2. Assessing eGRN enrichment in cell
## score the eRegulons in the cells to assess how enriched target regions and target genes are enriched in each cell

## Format eRegulons
from scenicplus.eregulon_enrichment import *
get_eRegulons_as_signatures(scplus_obj, eRegulon_metadata_key='eRegulon_metadata', key_added='eRegulon_signatures')

## Score chromatin layer
###  Region based raking
from scenicplus.cistromes import *
import time
start_time = time.time()
region_ranking = make_rankings(scplus_obj, target='region')

### Score region regulons
score_eRegulons(scplus_obj,
                ranking = region_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC', 
                enrichment_type= 'region',
                auc_threshold = 0.05,
                normalize = False,
                n_cpu = 1)
time = time.time()-start_time
print(time/60)

## Score transcriptome layer
### Gene based raking
from scenicplus.cistromes import *
import time
start_time = time.time()
gene_ranking = make_rankings(scplus_obj, target='gene')

### Score gene regulons
score_eRegulons(scplus_obj,
                gene_ranking,
                eRegulon_signatures_key = 'eRegulon_signatures',
                key_added = 'eRegulon_AUC', 
                enrichment_type = 'gene',
                auc_threshold = 0.05,
                normalize= False,
                n_cpu = 2)
time = time.time()-start_time
print(time/60)

# 3. Assessing TF-eGRN relationships

## Generate pseudobulks
import time
start_time = time.time()
generate_pseudobulks(scplus_obj, 
                         variable = 'ACC_cell_type',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Gene_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)

generate_pseudobulks(scplus_obj, 
                         variable = 'ACC_cell_type',
                         auc_key = 'eRegulon_AUC',
                         signature_key = 'Region_based',
                         nr_cells = 5,
                         nr_pseudobulks = 100,
                         seed=555)

time = time.time()-start_time
print(time/60)

## Correlation between TF and eRegulons
import time
start_time = time.time()
TF_cistrome_correlation(scplus_obj,
                        variable = 'ACC_cell_type', 
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Gene_based',
                        out_key = 'ACC_cell_type_eGRN_gene_based')
TF_cistrome_correlation(scplus_obj,
                        variable = 'ACC_cell_type', 
                        auc_key = 'eRegulon_AUC',
                        signature_key = 'Region_based',
                        out_key = 'ACC_cell_type_eGRN_region_based')
time = time.time()-start_time
print(time/60)

categories = sorted(set(scplus_obj.metadata_cell['ACC_cell_type']))

# 4. Identification of high quality regulons

df1 = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
df2 = scplus_obj.uns['eRegulon_AUC']['Region_based'].copy()
df1.columns = [x.split('_(')[0] for x in df1.columns]
df2.columns = [x.split('_(')[0] for x in df2.columns]
correlations = df1.corrwith(df2, axis = 0)
correlations = correlations[abs(correlations) > 0.6]

## Keep only R2G +
keep = [x for x in correlations.index if '+_+' in x] + [x for x in correlations.index if '-_+' in x] 

## Keep extended if not direct
extended = [x for x in keep if 'extended' in x]
direct = [x for x in keep if not 'extended' in x]
keep_extended = [x for x in extended if not x.replace('extended_', '') in direct]
keep = direct + keep_extended

## Keep regulons with more than 10 genes
keep_gene = [x for x in scplus_obj.uns['eRegulon_AUC']['Gene_based'].columns if x.split('_(')[0] in keep]
keep_gene = [x for x in keep_gene if (int(x.split('_(')[1].replace('g)', '')) > 10)]
keep_all = [x.split('_(')[0] for x in keep_gene]
keep_region = [x for x in scplus_obj.uns['eRegulon_AUC']['Region_based'].columns if x.split('_(')[0] in keep]
scplus_obj.uns['selected_eRegulons'] = {}
scplus_obj.uns['selected_eRegulons']['Gene_based'] = keep_gene
scplus_obj.uns['selected_eRegulons']['Region_based'] = keep_region

keep_gene # will return the list of high-quality regulon
