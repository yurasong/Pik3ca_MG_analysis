# Project directory
projDir = "./"

# Output directory
outDir = projDir + 'output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

# Temp dir
tmpDir = 'scratch/ysong/'
if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)

# 1. Creating a cisTopic object

from pycisTopic.cistopic_class import *

count_matrix=pd.read_csv(projDir+'count_table.tsv', sep='\t',index_col=0)
path_to_blacklist='input_files/mm10-blacklist.v2.bed' 

## Create cisTopic object
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix, path_to_blacklist=path_to_blacklist)

# Add cell information
# Adding cell information
cell_data =  pd.read_csv(projDir+'metadata.tsv', sep='\t',index_col=0)
cell_data

cistopic_obj.add_cell_data(cell_data)

# 2. Run model
## Comment: DO NOT run on the local workstation, it takes huge capacity of memory. 

path_to_mallet_binary="/home/audrey/Mallet/bin/mallet"
os.environ['MALLET_MEMORY'] = '60G' #set-up maximum memory for mallet -> if not needed, you can skip.

## Run models
models=run_cgs_models_mallet(path_to_mallet_binary,
                    cistopic_obj,
                    n_topics=[2,5,10,15,20,25,30,35,40],
                    n_cpu=8,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    tmp_path=tmpDir, #Use SCRATCH if many models or big data set
                    save_path=None)

## Save
import pickle
with open(outDir+'Mallet_models_500.pkl', 'wb') as f:
  pickle.dump(models, f)

# 3. Model selection

model=evaluate_models(models,
                     select_model=None,
                     return_model=True,
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= outDir + 'model_selection.pdf')

## Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)

## Save
with open(outDir + 'K8Pik_cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)

# 4. Clustering and visualisation

from pycisTopic.clust_vis import *
find_clusters(cistopic_obj,
                 target  = 'cell',
                 k = 10,
                 res = [0.8],
                 prefix = 'pycisTopic_',
                 scale = True,
                 split_pattern = '-')

## Running UMAP and tSNE on cistopic object 

run_umap(cistopic_obj,
                 target  = 'cell', scale=True)

run_tsne(cistopic_obj,
                 target  = 'cell', scale=True)

#os.mkdir(outDir+'/visualization')
plot_metadata(cistopic_obj,
                 reduction_name='UMAP',
                 variables=['cell_type', 'pycisTopic_leiden_10_0.8'], 
                 target='cell', num_columns=2,
                 text_size=10,
                 dot_size=5,
                 figsize=(15,5),
                 save= outDir + 'dimensionality_reduction_label.pdf')

## Plot topic-contributions

plot_topic(cistopic_obj,
            reduction_name = 'UMAP',
            target = 'cell',
            num_columns=5,
            save= outDir + 'dimensionality_reduction_topic_contr.pdf')

## Save
with open(outDir + 'K8Pik_cisTopicObject.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)

# 5. Topic binarization & qc

from pycisTopic.topic_binarization import *
region_bin_topics = binarize_topics(cistopic_obj, method='otsu', ntop=3000, plot=True, num_columns=5, save= outDir + 'otsu.pdf')

binarized_cell_topic = binarize_topics(cistopic_obj, target='cell', method='li', plot=True, num_columns=5, nbins=60)

from pycisTopic.topic_qc import *
topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

## Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
fig.savefig(outDir + 'Topic_qc.pdf', bbox_inches='tight')
plt.show()

## Annotation and matrix for each topic
topic_annot = topic_annotation(cistopic_obj, annot_var='cell_type', binarized_cell_topic=binarized_cell_topic, general_topic_thr = 0.2)
topic_qc_metrics = pd.concat([topic_annot[['cell_type', 'Ratio_cells_in_topic', 'Ratio_group_in_population']], topic_qc_metrics], axis=1)

## Save
with open(outDir + 'Topic_qc_metrics_annot.pkl', 'wb') as f:
  pickle.dump(topic_qc_metrics, f)
with open(outDir + 'binarized_cell_topic.pkl', 'wb') as f:
  pickle.dump(binarized_cell_topic, f)
with open(outDir + 'binarized_topic_region.pkl', 'wb') as f:
  pickle.dump(region_bin_topics, f)

# 6. Differentially Accessible Regions (DARs)

## Impute region accessibility exploiting the cell-topic and topic-region probabilities
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6) # shrink very low probability to 0

## log-normalise the impute data
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

## identify differentially accessible regions between groups.
## perform a Wilcoxon rank-sum test between each group in the specified variable and the rest
markers_dict= find_diff_features(cistopic_obj,
                      imputed_acc_obj,
                      variable='cell_type',
                      var_features=variable_regions,
                      contrasts=None,
                      adjpval_thr=0.05,
                      log2fc_thr=np.log2(1.5),
                      n_cpu=5,
                      split_pattern = '-')

# Save
with open(outDir + 'DARs/Imputed_accessibility.pkl', 'wb') as f:
  pickle.dump(imputed_acc_obj, f)
with open(outDir + 'DARs/DARs.pkl', 'wb') as f:
  pickle.dump(markers_dict, f)

# 7. Prepare cisTarget run

from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type', var_features=variable_regions, split_pattern = '-')


if not os.path.exists(os.path.join(projDir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(projDir, 'scATAC/candidate_enhancers'))

import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(projDir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(projDir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(projDir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))   