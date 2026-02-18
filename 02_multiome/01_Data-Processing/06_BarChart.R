#####################################################################################
# Bar plot of epithelial cell-type composition (K8Pik vs K8Pik_Klf5KO) – Fig. 5i
#
# Description:
#   Generates a grouped bar chart comparing epithelial cell-type proportions
#   between K8Pik and K8Pik_Klf5KO conditions. The script filters the input
#   dataset to selected samples and visualizes relative counts (ratio %)
#   per cell type using a publication-style ggplot2 layout.
#
# Workflow overview:
#     1. Load compositional data table
#     2. Subset to selected experimental conditions (K8Pik, K8Pik_Klf5KO)
#     3. Plot grouped bar chart with manual colour scheme
#
# Inputs:
#   - 06_barchart.csv
#       Required columns:
#         • CellType : annotated epithelial states
#         • count    : ratio (%) or normalized counts
#         • Sample   : experimental condition label
#
# Outputs:
#   - Grouped bar plot used for Fig. 5i
# Dependencies:
#   ggplot2, tidyverse
#
# Notes:
#   - `geom_col()` is used because values are precomputed.
#   - Axis text rotation improves readability for long cell-type labels.
#####################################################################################


# Library

library(ggplot2)
library(tidyverse)

# Data

data <- read.delim("06_barchart.csv", h=T, sep=",")

data_sel <- data %>% filter(Sample %in% c("K8Pik", "K8Pik_Klf5KO"))

# BarPlot (Fig. 5i)

ggplot(data_sel, aes(x = CellType, y = count, fill = Sample)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "black", linewidth = 0.2) +
  labs(x = NULL, y = "Ratio (%)", fill = "Condition") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 13)
  ) +
  scale_fill_manual(values=c("#00cf30", "#000fcc"))