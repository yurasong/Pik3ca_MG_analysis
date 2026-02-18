#####################################################################################
# Visualisation of prioritized ligand–receptor interactions: Extended Data Fig. 5n
#
# Description:
#   Generates a bubble plot summarising prioritized ligand–receptor (L–R)
#   interactions inferred from NicheNet analysis. Point size reflects
#   ligand activity (z-score, clipped for visualization) and colour encodes
#   the integrated prioritization score. Interactions are displayed by
#   signaling direction (CAF→HY vs HY→CAF).
#
# Workflow overview:
#     1. Load precomputed prioritisation results from Excel file
#     2. Clip activity z-score to improve dynamic range visualization
#     3. Generate interaction labels (ligand_receptor format)
#     4. Plot directional L–R activity using ggplot2 bubble plot
#
# Inputs:
#   - Input_plot.xlsx
#       Table containing ligand, receptor, activity_zscore,
#       prioritization_score, and Direction annotations
#
# Outputs:
#   - Bubble plot visualizing directional ligand–receptor interactions
#       Size  : L–R activity (activity_zscore_clipped)
#       Colour: prioritization score
#
#
# Dependencies:
#   ggplot2, xlsx
#
# Notes:
#   - Designed for CAF–HY signaling visualization in ERPik stromal analysis
#   - Uses viridis "plasma" colour scale for publication-ready aesthetics
#   - Axis titles intentionally removed for figure-panel integration
#####################################################################################


# Preparation
## Library

library(xlsx)
library(ggplot2)

## Data
df <- read.xlsx("Input_plot.xlsx", 1, h=T)

df$activity_zscore_clipped <- pmin(df$activity_zscore, 4)
df$interaction <- paste(df$ligand, "_", df$receptor, sep="")

# Plotting

ggplot(df, aes(x = interaction, y = Direction, size = activity_zscore_clipped,
               colour = prioritization_score)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = "L-R Activity", range = c(1, 6)) +
  scale_colour_viridis_c(name = "Prioritization score", option = "plasma", direction = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank()) +
  scale_y_discrete(limits=c("HY_CAF", "CAF_HY")) +
  scale_x_discrete(limits=c("TGFB2_TGFBR2", "TGFB2_TGFBR3", "TGFB2_TGFBR1", "TGFB2_ACVR1", "COL3A1_ITGB1",
                            "TGFB1_ACVR1", "TGFB1_ITGAV", "TGFB1_TGFBR3", "TGFB1_TGFBR2", "TGFB1_TGFBR1",
                            "COL15A1_ITGB1", "COL15A1_ITGA2", "COL3A1_ITGA2",
                            "TGFB1_ITGB6", "FGF7_FGFR2", "BMP2_BMPR1A", "BMP2_BMPR1B", 
                            "BMP2_ACVR2A", "BMP2_BMPR2", 
                            "FGF1_ITGB3", "FGF1_FGFR1", "FGF1_FGFR2", "FGF1_ITGAV"))
