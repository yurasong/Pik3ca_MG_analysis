#####################################################################################
# NicheNet ligand prioritisation analysis: CAF → HY signalling
#
# Description:
#   Performs ligand–receptor and ligand–target prioritisation using NicheNet
#   on a human-ortholog-converted Seurat object. CellPhoneDB v5 interactions
#   are reformatted and integrated with NicheNet v2 networks to infer
#   fibroblast (CAF)-to-hybrid epithelial (HY_BC/ER+) signaling programs.
#   The pipeline predicts active ligands, downstream targets, and prioritized
#   ligand–receptor interactions, followed by visualization and export.
#
# Workflow overview:
#     1. Load Seurat object with human gene symbols
#     2. Import and reformat CellPhoneDB v5 interaction database
#     3. Load NicheNet ligand–target and signaling networks
#     4. Define sender (Fibroblast) and receiver (HY_BC/ER+) populations
#     5. Identify expressed genes, receptors, and potential ligands
#     6. Compute ligand activity scores (AUPR-based prioritisation)
#     7. Extract ligand–target and ligand–receptor networks
#     8. Visualise heatmaps for activity, targets, and L–R interactions
#     9. Generate prioritisation tables integrating DE and expression metrics
#
# Inputs:
#   - ERPik_Stroma_human_converted.rds
#       Human-ortholog Seurat object derived from murine ERPik stroma dataset
#   - CellPhoneDB v5 interaction_input.csv (GitHub)
#   - NicheNet v2 network files (Zenodo):
#         ligand_target_matrix_nsga2r_final.rds
#         lr_network_human_21122021.rds
#         signaling_network_human_21122021.rds
#         gr_network_human_21122021.rds
#
# Outputs:
#   Directory: ./CAF-HY
#     - ligand_activities.xlsx          : ranked ligand activity scores
#     - ligand_prioritized.xlsx         : integrated ligand prioritisation table
#     - Heatmaps for ligand activity, ligand–target links, and L–R interactions
#ㅌㄴ
# Dependencies:
#   Seurat, nichenetr, tidyverse, dplyr, patchwork, xlsx,
#   httr, readr, tidyr, stringr, repr
#
# Notes:
#   - Gene symbols are assumed to be human orthologs prior to analysis
#   - CellPhoneDB interactions are reformatted to match NicheNet structure
#   - Failed ligands during activity prediction are automatically removed
#   - Custom visualisation relies on script/mushroom_plot.R
#####################################################################################

# Preparation
## Library

library(httr)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(repr)

library(nichenetr)
library(Seurat) 
library(tidyverse) 
library(patchwork)
library("xlsx")

## Data
seuset <- readRDS("ERPik_Stroma_human_converted.rds") # Converted to Human genes

# Configuration

options(timeout = 3600)
options(repr.plot.width=15, repr.plot.height=10)

dir.create("CAF-HY")
saveLoc <- "./CAF-HY"

source("script/mushroom_plot.R")

## Database
### CPDB v5 database

url <- "https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv" 
resource <- read_csv(url)

resource <- resource %>%
  mutate(interactors = str_replace_all(interactors, "\\+", "_")) # Replace '+' with '_'

resource <- resource %>%
  mutate(interactors = if_else(str_count(interactors, "-") == 2,
                               sub("-", "&", interactors),
                               interactors)) # If 'interactors' contains two '-', replace the first one with '&'

interactors_split <- str_split_fixed(resource$interactors, "-", 2) # Split 'interactors' by '-' and expand into two columns

ligand <- str_replace_all(interactors_split[,1], "&", "-") # Replace '&' back to '-' in the first column
receptor <- interactors_split[,2]

resource <- tibble(from = ligand, to = receptor, database = "CellphoneDBV5", source = "CellphoneDBV5")

### Ligand-target matrix

zenodo_path <- "https://zenodo.org/record/7074291/files/"
ligand_target_matrix <- readRDS(url(paste0(zenodo_path, "ligand_target_matrix_nsga2r_final.rds")))

lr_network <- resource %>%  separate_rows(to, sep = "_") %>%
  separate_rows(from, sep = "_") %>%
  distinct()

weighted_networks <- readRDS(url(paste0(zenodo_path, "weighted_networks_nsga2r_final.rds"))) 

length(unique(lr_network$to)) # Should return 669

lr_network <- lr_network[lr_network$from %in% (colnames(ligand_target_matrix)),]
length(unique(lr_network$to)) # Should return 551

### Human network data: Download related files from NicheNet-v2

zenodo_path <- "https://zenodo.org/record/7074291/files/"
lr_network_human <- readRDS(url(paste0(zenodo_path, "lr_network_human_21122021.rds")))
sig_network_human <- readRDS(url(paste0(zenodo_path, "signaling_network_human_21122021.rds")))
gr_network_human <- readRDS(url(paste0(zenodo_path, "gr_network_human_21122021.rds")))

weighted_networks <- construct_weighted_networks(
  lr_network = lr_network,
  sig_network = sig_network_human,
  gr_network = gr_network_human,
  source_weights_df = source_weights)

# NicheNet analysis: CAF to HY

receiver <- "HY_BC/ER+" 
expressed_genes_receiver <- get_expressed_genes(receiver, seuset,  pct = 0.1)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver) 

potential_ligands <- lr_network[lr_network$to %in% expressed_receptors, ] 
potential_ligands <- unique(potential_ligands$from) 

sender_celltypes <- c("Fibroblast") 
list_expressed_genes_sender <- lapply(sender_celltypes, function(celltype) {
  get_expressed_genes(celltype, seuset, pct = 0.1)
}) 
expressed_genes_sender <- unique(unlist(list_expressed_genes_sender)) 
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

DE_table_receiver <- FindMarkers(object = seuset,  
                                 ident.1 = "HY_BC/ER+",
                                 ident.2 = "Fibroblast",
                                 group.by = "cell_type",
                                 min.pct = 0.1) 

geneset_oi <- DE_table_receiver[DE_table_receiver$p_val_adj <= 0.01 &
                                  (DE_table_receiver$avg_log2FC) >= 0.5, ] 
geneset_oi <- rownames(geneset_oi)[rownames(geneset_oi) %in% rownames(ligand_target_matrix)] 

background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in% rownames(ligand_target_matrix)] 

## To filter out which is not in ligand

fails <- character(0)

for (lig in potential_ligands) {
  ok <- tryCatch({
    predict_ligand_activities(
      geneset = geneset_oi,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = lig
    )
    TRUE
  }, error = function(e) {
    message("FAILED ligand: ", lig, " | ", conditionMessage(e))
    FALSE
  })
  
  if (!ok) fails <- c(fails, lig)
}

potential_ligands <- setdiff(potential_ligands, fails)

## L-R Activity calculation and export

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities[order(ligand_activities$aupr_corrected, decreasing = TRUE), ] 

ligand_activities_all <- ligand_activities 
ligand_activities <- ligand_activities[ligand_activities$test_ligand %in% potential_ligands_focused, ] 

write.xlsx(ligand_activities, paste0(saveLoc, "/ligand_activities.xlsx"), col.names = TRUE, row.names = TRUE, append = FALSE)

## Best L-R pair

best_upstream_ligands <- top_n(ligand_activities, 50, aupr_corrected)$test_ligand 

active_ligand_target_links_df <- lapply(best_upstream_ligands,
                                        get_weighted_ligand_target_links, 
                                        geneset = geneset_oi, 
                                        ligand_target_matrix = ligand_target_matrix, 
                                        n = 200) 

active_ligand_target_links_df <- drop_na(bind_rows(active_ligand_target_links_df)) 

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands = best_upstream_ligands,
  expressed_receptors = expressed_receptors,
  lr_network = lr_network,
  weighted_networks_lr_sig = weighted_networks$lr_sig
)

## Visualisation of the results: Used for sanity check of the analysis

ligand_aupr_matrix <- column_to_rownames(ligand_activities, "test_ligand") 
ligand_aupr_matrix <- ligand_aupr_matrix[rev(best_upstream_ligands), "aupr_corrected", drop=FALSE] 
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p2= (make_heatmap_ggplot(vis_ligand_aupr,
                         y_name = "Prioritized ligands", x_name = "Ligand activity",
                         legend_title = "AUPR", color = "darkorange") + 
       theme(axis.text.x.top = element_blank(), text = element_text(size = 13)))
p2

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25) 

order_ligands <- rev(intersect(best_upstream_ligands, colnames(active_ligand_target_links))) 
order_targets <- intersect(unique(active_ligand_target_links_df$target), rownames(active_ligand_target_links)) 

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

(make_heatmap_ggplot(vis_ligand_target,
                     y_name = "Prioritized ligands", x_name = "Predicted target genes",
                     color = "purple", legend_title = "Regulatory potential") + 
    scale_fill_gradient2(low = "whitesmoke",  high = "purple")) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df, best_upstream_ligands,
  order_hclust = "receptors") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Prioritized ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential") 
  + theme(text = element_text(size = 15)))

## Prioritisation

lr_network_filtered <- filter(lr_network, from %in% potential_ligands_focused &
                                to %in% expressed_receptors)[, c("from", "to")]

info_tables <- generate_info_tables( 
  seuset, 
  celltype_colname = "cell_type", 
  senders_oi = "Fibroblast", 
  receivers_oi = "HY_BC/ER+", 
  lr_network = lr_network_filtered, 
  scenario = "one_condition"
) 

processed_DE_table <- info_tables$sender_receiver_de  
processed_expr_table <- info_tables$sender_receiver_info  
processed_condition_markers <- info_tables$lr_condition_de

head(processed_DE_table,5)

prioritized_table <- generate_prioritization_tables(
  sender_receiver_info = processed_expr_table,
  sender_receiver_de = processed_DE_table,
  ligand_activities = ligand_activities,
  lr_condition_de = processed_condition_markers,
  scenario = "one_condition") 

prioritized_table$sender <- factor(prioritized_table$sender, levels = sender_celltypes)

prioritizing_weights <- c("de_ligand" = 1,
                          "de_receptor" = 1,
                          "activity_scaled" = 1,
                          "exprs_ligand" = 1,
                          "exprs_receptor" = 1,
                          "ligand_condition_specificity" = 1,
                          "receptor_condition_specificity" = 1) 

prioritized_table$sender <- factor(prioritized_table$sender, levels = sender_celltypes)

write.xlsx(prioritized_table, paste0(saveLoc, "/ligand_prioritized.xlsx"), col.names = TRUE, row.names = TRUE, append = FALSE)
