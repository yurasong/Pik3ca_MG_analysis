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
