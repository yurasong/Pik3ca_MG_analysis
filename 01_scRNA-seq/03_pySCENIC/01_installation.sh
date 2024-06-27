#!/bin/bash

############################################################
# This run has been proceeded on Ubuntu 18.04.
# All can be done on terminal. This is shell based.
############################################################

# It is always recommended to install the pyscenic with anaconda, since one of the dependencies always makes the crash.
# Please use python version at least 3.7, the previous versions may not install the dependencies of pyscenic properly.

conda create -y -n pyscenic python=3.7
conda activate pyscenic
conda install -y numpy
conda install -y -c anaconda cytoolz

pip install pyscenic

########### DATABASE SETTING ##############

# For getting transcription factor information which is compatible with the pyscenic:
## When this installation works, you will have pySCENIC directory under the home directory. You could find files under the resources directory.

git clone https://github.com/aertslab/pySCENIC.git

# feature file download: Also you could download manually by connecting to https://resources.aertslab.org/cistarget/databases/

## Human database:
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feature

wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feature

## Mouse database:

wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feature

wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feature

## TF annotation:

wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
