#####################################################################################
# Bar plot visualisation of cell-type proportions across conditions (Fig. 5g)
#
# Description:
#   Generates a grouped bar chart showing the percentage distribution of
#   epithelial cell states across experimental conditions. Cell types are
#   displayed in a predefined biological order, and bars are colour-coded
#   by sample/condition. This plot is typically used for figure panels
#   summarising compositional changes between ER Pik–derived populations.
#
# Workflow overview:
#     1. Load tabular data containing cell_type, percentage, and Sample columns
#     2. Set factor levels to enforce biologically meaningful ordering
#     3. Plot grouped bar chart using ggplot2 with manual colour scheme
#
# Inputs:
#   - barchart.csv
#       Columns required:
#         • cell_type   : annotated epithelial states
#         • percentage  : ratio (%) per state
#         • Sample      : experimental condition label
#
# Outputs:
#   - Grouped bar plot showing ratio (%) per cell type and condition
#
#
# Dependencies:
#   ggplot2, tidyverse
#
# Notes:
#   - `geom_col()` is used since values are pre-computed percentages.
#   - Manual factor ordering ensures consistent figure layout across panels.
#####################################################################################



# Library

library(ggplot2)
library(tidyverse)

# Data

data <- read.delim("barchart.csv", h=T, sep=",")

data$cell_type <- factor(
  data$cell_type,
  levels = c(
    "ER+ LCs",
    "Early HY",
    "Late HY",
    "Myo",
    "HY ER+/ER-"
  ))

# BarPlot (Fig. 5g)

ggplot(data, aes(x = cell_type, y = percentage, fill = Sample)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "black", linewidth = 0.2) +
  labs(x = NULL, y = "Ratio (%)", fill = "Condition") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 13)
  ) +
  scale_fill_manual(values=c("#00cf30", "#000fcc"))