#####################################################################################
# File: 06_PieChart.R
#
# Description:
#   Reads cell‐type count ratios for CTL, K8Pik, and K8Pik_Klf5KO samples,
#   and generates pie charts for each sample.
#
# Inputs:
#   - input_data/Chart_ratio.csv : CSV file with columns Sample, CellType, count
#
# Outputs: Corresponding to Fig. 6d
#   - PieChart_CTL.pdf
#   - PieChart_K8Pik.pdf
#   - PieChart_K8Pik_Klf5KO.pdf
#
# Annotation was done manually after generating piechart.
#
# Dependencies:
#   dplyr, data.table, tidyverse, RColorBrewer, ggplot2
#####################################################################################

# Library
library(dplyr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# Data
ratio <- read.delim("input_data/Chart_ratio.csv", header = TRUE, sep = ",")

ctl <- ratio %>% filter(Sample == "CTL")
pik <- ratio %>% filter(Sample == "K8Pik")
ko  <- ratio %>% filter(Sample == "K8Pik_Klf5KO")

# Plotting
color_palette <- c(
  "#CD9600", "#F8766D", "#7CAE00", "#00BE67",
  "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"
)

plot_pie <- function(data, filename) {
  p <- ggplot(data, aes(x = "", y = cou


# Library
library(dplyr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# Data
ratio <- read.delim("input_data/Chart_ratio.csv", h=T, sep=",")

ctl <- data %>% filter(Sample == "CTL")
pik <- data %>% filter(Sample =="K8Pik")
ko <- data %>% filter(Sample =="K8Pik_Klf5KO")

ggplot(ctl, aes(x="", y=count, fill=CellType)) +
  geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0) + 
  theme_void() + 
  scale_fill_manual(values=c("#CD9600", "#F8766D", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"))

ggsave("PieChart_CTL.pdf", height=10, width=10)

ggplot(pik, aes(x="", y=count, fill=CellType)) +
  geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0) + 
  theme_void() + 
  scale_fill_manual(values=c("#CD9600", "#F8766D", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"))

ggsave("PieChart_K8Pik.pdf", height=10, width=10)

ggplot(ko, aes(x="", y=count, fill=CellType)) +
  geom_bar(stat="identity", width=1, color="white") + 
  coord_polar("y", start=0) + 
  theme_void() + 
  scale_fill_manual(values=c("#CD9600", "#F8766D", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"))

ggsave("PieChart_K8Pik_Klf5KO.pdf", height=10, width=10)