#!/bin/bash

############################################################
# This run has been proceeded on Ubuntu 18.04.
# All can be done on terminal. This is shell based.
############################################################

# Create bamfile which only includes the cells passed QC
# This will return the bam file per cell type.
sinto filterbarcodes -b K8Pik_atac_possorted_bam.bam -c 00_input/metadata_K8Pik_noFibro_multiome.txt
sinto filterbarcodes -b K8Pik_atac_possorted_bam.bam -c 00_input/metadata_K8Pik_Klf5KO_noFibro_multiome.txt

# Filter out unmapped contigs
## Here, I am taking the Immature_BC.bam of Klf5KO sample as an example. The procedure is exactly same for the other bam files.

samtools view -bh Immature_BC.bam chr1 > chr1.bam 
samtools view -bh Immature_BC.bam chr2 > chr2.bam 
samtools view -bh Immature_BC.bam chr3 > chr3.bam 
samtools view -bh Immature_BC.bam chr4 > chr4.bam 
samtools view -bh Immature_BC.bam chr5 > chr5.bam 
samtools view -bh Immature_BC.bam chr6 > chr6.bam 
samtools view -bh Immature_BC.bam chr7 > chr7.bam 
samtools view -bh Immature_BC.bam chr8 > chr8.bam 
samtools view -bh Immature_BC.bam chr9 > chr9.bam 
samtools view -bh Immature_BC.bam chr10 > chr10.bam 
samtools view -bh Immature_BC.bam chr11 > chr11.bam 
samtools view -bh Immature_BC.bam chr12 > chr12.bam 
samtools view -bh Immature_BC.bam chr13 > chr13.bam 
samtools view -bh Immature_BC.bam chr14 > chr14.bam 
samtools view -bh Immature_BC.bam chr15 > chr15.bam 
samtools view -bh Immature_BC.bam chr16 > chr16.bam 
samtools view -bh Immature_BC.bam chr17 > chr17.bam 
samtools view -bh Immature_BC.bam chr18 > chr18.bam 
samtools view -bh Immature_BC.bam chr19 > chr19.bam 
samtools view -bh Immature_BC.bam chrX > chrX.bam 
samtools view -bh Immature_BC.bam chrY > chrY.bam

samtools merge -O BAM Immature_BC_noContig.bam chr*.bam
samtools sort -O BAM -T tempsort -o Immature_BC_noContig_sorted.bam Immature_BC_noContig.bam
samtools index -b Immature_BC_noContig_sorted.bam Immature_BC_noContig_sorted.bam.bai

samtools view -b -F 4 Immature_BC_noContig_sorted.bam > mapped.bam_file.bam
samtools view -bf 1 mapped.bam_file.bam > paired_end_mapped.bam

