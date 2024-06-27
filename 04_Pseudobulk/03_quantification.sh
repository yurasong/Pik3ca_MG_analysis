#!/bin/bash

############################################################
# This run has been proceeded on Ubuntu 18.04.
# All can be done on terminal. This is shell based.
# The code has been run on Immature_BC.bam of Klf5KO sample as an example.
# Proceduce is exactly same for the other bam files.
############################################################

# Generate merged peak profile
cat blacklist_removed_* >> total.bed # Here use all the blacklisted removed bed file from each bamfile.
bedtools sort -i total.bed > total_sorted.bed
bedtools merge -i total_sorted.bed > merge.bed

awk '{FS="\t";OFS="\t";$4="merged"NR;print}' merge.bed > merge1.bed
awk '{FS="\t";OFS="\t";$6=".";$7="+";$8=".";$9="gene_id \""$4"\"";$4=$2;$5=$3;$2="merged";$3="exon";print}' merge1.bed > merge.gff

# Quantification
htseq-count -f bam -r pos -m Immature_BC/Klf5KO/paired_end_mapped.bam merge.gff > count_Klf5KO_Immature_BC.txt