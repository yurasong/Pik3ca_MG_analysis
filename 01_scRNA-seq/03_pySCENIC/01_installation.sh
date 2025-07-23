#!/bin/bash

####################################################################################################################################
# 01_installation.sh
# Installs pyscenic via conda + pip to avoid crash between dependencies
# Python version >= 3.7 required
# downloads the required cistarget databases
####################################################################################################################################

conda create -y -n pyscenic python=3.7
conda activate pyscenic
conda install -y numpy
conda install -y -c anaconda cytoolz

pip install pyscenic

########### DATABASE SETTING ##############

# For getting transcription factor information which is compatible with the pyscenic:
## When this installation works, you will have pySCENIC directory under the home directory. You could find files under the resources directory.

git clone https://github.com/aertslab/pySCENIC.git

########### If prefer manual download ##############

# feature file available via https://resources.aertslab.org/cistarget/databases/

## Human database:
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feature
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feature

## Mouse database:

wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feature
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feature

## TF annotation:

wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
