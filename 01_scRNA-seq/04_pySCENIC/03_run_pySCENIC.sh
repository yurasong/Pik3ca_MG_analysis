#!/bin/bash

########################################################################################################################
# 03_run_scenic.sh
# GRN → ctx → aucell workflows in a loop for N replicates.
# Run on pyscenic conda environment: conda activate pyscenic before running this code!
########################################################################################################################

EXPR_MTX="exprMat.tsv" # Expression matrix from scRNA-seq data
TF_LIST="/home/audrey/pySCENIC/resources/mm_mgi_tfs.txt"
DB_GLOB="/home/audrey/pySCENIC/cisTarget_databases/mm10/mm10*"
ANNOT="/home/audrey/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
N_WORKERS=10 # For parallel
N_RUNS=10 # Number of jobs to repeat, to correct stochastic effect

# ——— Main loop —————————————————————————————
for (( i=0; i< N_RUNS; i++ )); do
  echo "=== Starting run #${i} ==="

  # 1) Inference of co‐expression modules (GRN)
  OUT_ADJ="adj${i}.csv"
  echo "  → pyscenic grn (output: ${OUT_ADJ})"
  pyscenic grn \
    --num_workers "${N_WORKERS}" \
    -o "${OUT_ADJ}" \
    "${EXPR_MTX}" "${TF_LIST}"

  # 2) Motif enrichment (ctx)
  OUT_REG="regulons${i}.csv"
  echo "  → pyscenic ctx (output: ${OUT_REG})"
  pyscenic ctx \
    "${OUT_ADJ}" "${DB_GLOB}" \
    --annotations_fname "${ANNOT}" \
    --expression_mtx_fname "${EXPR_MTX}" \
    --mode "dask_multiprocessing" \
    --output "${OUT_REG}" \
    --num_workers "${N_WORKERS}"

  # 3) AUCell scoring
  OUT_AUC="auc_mtx${i}.csv"
  echo "  → pyscenic aucell (output: ${OUT_AUC})"
  pyscenic aucell \
    "${EXPR_MTX}" "${OUT_REG}" \
    -o "${OUT_AUC}" \
    --num_workers "${N_WORKERS}"

  echo "=== Finished run #${i} ==="$'\n'
done

echo "All ${N_RUNS} runs complete."