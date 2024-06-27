#####################################################################################
# This code is applied on Integrated data.
#####################################################################################

# Library
library(dplyr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# Data
ratio <- read.delim("input_data/Chart_ratio.csv", h=T, sep=",")

# Plotting
ggplot(data=ratio, aes(x=CellType, y=count, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge())+
  theme_minimal() +
  scale_x_discrete(limits=c("LC_ER+", "HY_ER+_ER-", "LC_ER-", "HY_BC/ER-", "Immature_BC",
                            "Myoepith", "HY_BC/ER+", "Prolif")) +
  scale_fill_manual(values=c('#FA8072','#008B8B')) +
  geom_text(aes(label=count), vjust=0, color="black",
            position=position_dodge(0.9), size=3.5) +
  ylab("Percentage (%)") +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.title.x=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"))