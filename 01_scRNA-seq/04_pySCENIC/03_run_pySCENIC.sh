#!/bin/bash

############################################################
# This run has been proceeded on Ubuntu 18.04.
# All can be done on terminal. This is shell based.
# This script takes the expression matrix from seurat object as an input.
# Not going through the loom file creating and building annotation files. 
# Please follow 01_installation.sh file before starting SCENIC.
############################################################

conda activate pyscenic

pyscenic grn --num_workers 10 -o adj0.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj0.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons0.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons0.csv -o auc_mtx0.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj1.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj1.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons1.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons1.csv -o auc_mtx1.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj2.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj2.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons2.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons2.csv -o auc_mtx2.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj3.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj3.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons3.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons3.csv -o auc_mtx3.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj4.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj4.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons4.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons4.csv -o auc_mtx4.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj5.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj5.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons5.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons5.csv -o auc_mtx5.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj6.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj6.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons6.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons6.csv -o auc_mtx6.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj7.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj7.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons7.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons7.csv -o auc_mtx7.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj8.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj8.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons8.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons8.csv -o auc_mtx8.csv --num_workers 10

pyscenic grn --num_workers 10 -o adj9.csv exprMat.tsv /home/audrey/pySCENIC/resources/mm_mgi_tfs.txt 
pyscenic ctx adj9.csv /home/audrey/pySCENIC/cisTarget_databases/mm10/mm10* --annotations_fname /home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons9.csv --num_workers 10
pyscenic aucell exprMat.tsv regulons9.csv -o auc_mtx9.csv --num_workers 10
