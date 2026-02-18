#!/bin/bash

##########################################################################################
# File: 04_bedgraph_generation.sh
#
# Description:
#   Scales per-cell-type BAM coverage to bedGraph format using custom scale factors.
#   For each entry (cellDir, sample, scale):
#     - Reads the paired-end mapped BAM
#     - Generates a bedGraph with genome coverage scaled by the provided factor
#
# Note:
#   Scale factor is calculated as (normalized into 1 M reads) / (raw count) for each sample.
#
# Configuration:
#   - entries: array of "<cellDir> <sample> <scale>" triples
#
# Inputs:
#   - ${cellDir}/${sample}/paired_end_mapped.bam  : Paired-end, mapped BAM per cell type
#
# Outputs:
#   - scaled_<sample>_<cellDir>.bedgraph          : Scaled coverage bedGraph files
#
# Dependencies:
#   - bash
#   - bedtools (genomecov)
##########################################################################################

#----------------------------------------
# COVERAGE SCALING CONFIG
#----------------------------------------
# List of triples: "<cellDir> <sample> <scale>"
entries=(
  "HY_BC_ER+        K8Pik    11.71"
  "HY_BC_ER+        Klf5KO   33.40"
  "HY_ER+_ER-       K8Pik     4.29"
  "HY_ER+_ER-       Klf5KO   46.00"
  "Immature_BC      K8Pik     1.49"
  "Immature_BC      Klf5KO    2.97"
  "LC_ER+           K8Pik     1.78"
  "LC_ER+           Klf5KO    0.90"
  "LC_ER-           K8Pik     3.76"
  "LC_ER-           Klf5KO    1.00"
  "K8PIk_only/Myoepith    K8Pik    4.25"
  "K8PIk_only/HY_BC_ER-   K8Pik    3.57"
)

#----------------------------------------
# MAIN LOOP
#----------------------------------------
for entry in "${entries[@]}"; do
  read -r cellDir sample scale <<<"$entry"

  bam="${cellDir}/${sample}/paired_end_mapped.bam"
  # build a safe output name: replace any “/” in cellDir with “_”
  safeCell="${cellDir//\//_}"
  out="scaled_${sample}_${safeCell}.bedgraph"

  if [[ ! -f "$bam" ]]; then
    echo "Missing BAM: $bam – skipping"
    continue
  fi

  echo "▶ Scaling coverage for ${sample} / ${cellDir} (×${scale}) → ${out}"
  bedtools genomecov \
    -ibam "$bam" \
    -bg \
    -scale "$scale" \
    > "$out"
done

echo "All bedgraph scaling complete."