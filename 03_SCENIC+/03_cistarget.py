#####################################################################################
# File: 03_pycistarget.py
#
# Description:
#   1. Load binarized topic region sets and DARs (OTSU, top 3k, and DARs) from pickle files
#   2. Filter for standard chromosomes and convert region names to genomic coordinates
#   3. Organize regions into PyRanges objects for each topic/DAR
#   4. Run pycistarget motif enrichment on each region set using specified rankings/scores databases
#   5. Save output motif files under the `motifs/` directory
#
#
# Inputs (relative to projDir):
#   - scATAC/candidate_enhancers/region_bin_topics_otsu.pkl
#   - scATAC/candidate_enhancers/region_bin_topics_top3k.pkl
#   - scATAC/candidate_enhancers/markers_dict.pkl
#   - ctx_db_path (feather): mm10_screen_v10_clust.regions_vs_motifs.rankings.feather
#   - dem_db_path (feather): mm10_screen_v10_clust.regions_vs_motifs.scores.feather
#   - motif_annotation (tbl): motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
#
# Outputs (under projDir/motifs/):
#   - One motif enrichment result file per region set (e.g., topics_otsu_<topic>.feather, 
#     topics_top_3_<topic>.feather, DARs_<DAR>.feather)
#
# Dependencies:
#   os, pickle, pandas, pyranges, scenicplus, pycistarget.utils.region_names_to_coordinates,
#   scenicplus.wrappers.run_pycistarget
#
# Notes:
#   - Ensure the MALLET_MEMORY and other environment variables are set if required.
#   - Adjust n_cpu and _temp_dir parameters in run_pycistarget call as needed.
#####################################################################################


import sys
import os

# Project directory
projDir = "/globalscratch/ysong/"

# Output directory
outDir = projDir + 'output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

import pickle
region_bin_topics_otsu = pickle.load(open(os.path.join(projDir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(projDir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(projDir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')

## Run cistarget

db_fpath = "/home/ulb/iribhm/ysong/database_scenicplus"
rankings_db = os.path.join(db_fpath, 'mm10_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'mm10_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(db_fpath, 'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl')

if not os.path.exists(os.path.join(projDir, 'motifs')):
    os.makedirs(os.path.join(projDir, 'motifs'))

from scenicplus.wrappers.run_pycistarget import run_pycistarget

run_pycistarget(
    region_sets = region_sets,
    species = 'mus_musculus',
    save_path = os.path.join(projDir, 'motifs'),
    ctx_db_path = rankings_db,
    dem_db_path = scores_db,
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 4, _temp_dir = "/scratch/ysong", 
    annotation_version = 'v10nr_clust',
    )
