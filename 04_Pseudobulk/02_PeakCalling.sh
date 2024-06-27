#!/bin/bash

############################################################
# This run has been proceeded on Ubuntu 18.04.
# All can be done on terminal. This is shell based.
# The code has been run on Immature_BC.bam of Klf5KO sample as an example.
# Proceduce is exactly same for the other bam files.
############################################################

# Peak calling with MACS2
macs2 callpeak -t Immature_BC_noContig_sorted.bam --outdir Klf5KO_immature_BC -f BAMPE -g mm -n Klf5KO_immature_BC -q 0.01 --nomodel

## blacklisted region removal
cd Klf5KO_immature_BC
sed '/^$/d' *.xls > peaks.bed
sed '/#/d' peaks.bed > peaks1.bed
sed '1d' peaks1.bed > peaks.bed
awk '{ print $1"\t"$2"\t"$3"\t"$10"\t"$7}' peaks.bed > peaks1.bed
mv peaks1.bed peaks.bed

bedtools intersect -v -a peaks.bed -b 00_input/mm10-blacklist.v2.bed > blacklist_removed_Klf5KO_immature_BC.bed
rm peaks.bed






