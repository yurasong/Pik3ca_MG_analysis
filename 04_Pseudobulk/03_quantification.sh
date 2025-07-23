#!/bin/bash

##########################################################################################
# File: 03_quantification.sh
#
# Description:
#   1. Merge per-cell-type peak BED files into a unified peak set:
#       - Concatenate all blacklist-removed BEDs
#       - Sort and merge overlapping peaks
#       - Generate a GFF annotation file for merged peaks
#   2. Quantify read counts per merged peak for each cell type using htseq-count
#
# Configuration:
#   - annotation: path to merged GFF file ("merge.gff")
#   - metas: list of metadata files with sample and cell-type columns
#
# Inputs:
#   - {sample}/{cell}_paired.bam       : filtered, paired BAMs per cell type
#   - merge.gff                        : GFF of merged peaks
#   - 00_input/metadata_<sample>.txt   : metadata listing cell types per sample
#
# Outputs:
#   - merge.bed, merge.gff             : merged peak definitions
#   - {sample}/count_<sample>_<cell>.txt : peak count tables per cell type
#
# Dependencies:
#   - bash
#   - bedtools
#   - awk, sed, cat
#   - htseq-count (from HTSeq)
##########################################################################################

# Generate merged peak profile
cat blacklist_removed_* >> total.bed # Here use all the blacklisted removed bed file from each bamfile.
bedtools sort -i total.bed > total_sorted.bed
bedtools merge -i total_sorted.bed > merge.bed

awk '{FS="\t";OFS="\t";$4="merged"NR;print}' merge.bed > merge1.bed
awk '{FS="\t";OFS="\t";$6=".";$7="+";$8=".";$9="gene_id \""$4"\"";$4=$2;$5=$3;$2="merged";$3="exon";print}' merge1.bed > merge.gff

#----------------------------------------
# CONFIGURATION (add to top of your script)
#----------------------------------------
# Path to your merged GFF
annotation="merge.gff"

# reuse your existing metadata list
metas=(
  "00_input/metadata_K8Pik_noFibro_multiome.txt"
  "00_input/metadata_K8Pik_Klf5KO_noFibro_multiome.txt"
  "00_input/metadata_CTL_multiome.txt"
)

#----------------------------------------
# FUNCTION: run htseq-count for one cell-type
#   $1 = sample name
#   $2 = cell-type name
#----------------------------------------
quantify_bam() {
  local sample=$1
  local cell=$2
  local bam="${sample}/${cell}_paired.bam"
  local out="${sample}/count_${sample}_${cell}.txt"

  if [[ ! -f "$bam" ]]; then
    echo "Missing $bam, skipping"
    return
  fi

  echo "▶ [$sample/$cell] htseq-count → $out"
  htseq-count \
    -f bam \
    -r pos \
    -m intersection-nonempty \
    "$bam" \
    "$annotation" \
    > "$out"
  echo "✔ Wrote $out"
}

#----------------------------------------
# MAIN (append after your contig & peak steps)
#----------------------------------------
echo "=== Starting quantification with htseq-count ==="
for meta in "${metas[@]}"; do
  sample=$(basename "$meta")
  sample=${sample#metadata_}
  sample=${sample%_multiome.txt}

  echo "--- Sample: $sample ---"
  # read unique cell types from 2nd column
  while read -r cell; do
    quantify_bam "$sample" "$cell"
  done < <(cut -f2 "00_input/${meta##*/}" | sort -u)
done
echo "All quantification complete."
