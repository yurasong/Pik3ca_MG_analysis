#!/bin/bash

##########################################################################################
# File: 02_PeakCalling.sh
#
# Description:
#   For each sample and cell type:
#     1) Call peaks on contig-stripped, sorted BAMs using MACS2
#     2) Extract and format peak coordinates from MACS2 XLS output
#     3) Remove blacklist regions from the peak list using Bedtools
#
# Configuration (edit within script):
#   - metas array: paths to per-sample metadata files listing cell types
#   - blacklist: path to mm10 blacklist BED file
#
# Inputs (per sample directory):
#   - {cell}_noContig_sorted.bam      : filtered & sorted BAM for each cell type
#   - 00_input/metadata_{sample}_multiome.txt : metadata listing cell types
#   - mm10-blacklist.v2.bed            : genomic blacklist regions
#
# Outputs:
#   - {sample}_{cell}/                 : MACS2 output directory per cell type
#       * peaks.bed                    : formatted peak BED
#       * blacklist_removed_{name}.bed : final peak set without blacklist regions
#
# Dependencies:
#   - bash
#   - MACS2 (macs2 callpeak)
#   - Bedtools (bedtools intersect)
#   - coreutils (sed, awk, mkdir, cd)
##########################################################################################

#----------------------------------------
# CONFIGURATION
#----------------------------------------
# (reuse your existing metas array)

metas=(
  "00_input/metadata_K8Pik_noFibro_multiome.txt"
  "00_input/metadata_K8Pik_Klf5KO_noFibro_multiome.txt"
  "00_input/metadata_CTL_multiome.txt"
)

# path to mm10 blacklist
blacklist="00_input/mm10-blacklist.v2.bed"

#----------------------------------------
# FUNCTION: call MACS2 + clean peaks
#   $1 = sample name
#   $2 = cell type
#----------------------------------------
call_and_clean_peaks() {
  local sample=$1
  local cell=$2
  local bam="${cell}_noContig_sorted.bam"
  local name="${sample}_${cell}"
  local outdir="${name}"

  echo "▶ [${sample}] calling peaks on ${bam}"
  macs2 callpeak \
    -t "$bam" \
    --outdir "$outdir" \
    -f BAMPE \
    -g mm \
    -n "$name" \
    -q 0.01 \
    --nomodel

  echo "▶ [${name}] removing blacklist regions"
  (
    cd "$outdir" || exit
    # pull out the raw XLS peak list
    sed '/^$/d' *.xls       > peaks.bed
    sed '/#/d' peaks.bed     > peaks1.bed
    sed '1d' peaks1.bed      > peaks.bed
    awk '{ print $1"\t"$2"\t"$3"\t"$10"\t"$7 }' peaks.bed \
                             > peaks1.bed
    mv peaks1.bed peaks.bed

    # intersect‐v to drop blacklist
    bedtools intersect -v \
      -a peaks.bed \
      -b "../../$blacklist" \
      > "blacklist_removed_${name}.bed"

    rm peaks.bed
  )
  echo "Done peaks for ${name}"
}

#----------------------------------------
# MAIN: iterate samples → cell types → peaks
#----------------------------------------
for meta in "${metas[@]}"; do
  sample=$(basename "$meta")
  sample=${sample#metadata_}
  sample=${sample%_multiome.txt}

  echo "=== Peak-calling for sample: ${sample} ==="
  pushd "$sample" >/dev/null || { echo " No dir $sample, skipping"; continue; }

  # for each unique cell type in column 2 of the metadata:
  while read -r cell; do
    if [[ -f "${cell}_noContig_sorted.bam" ]]; then
      call_and_clean_peaks "$sample" "$cell"
    else
      echo " Missing ${cell}_noContig_sorted.bam, skipping"
    fi
  done < <(cut -f2 "../$meta" | sort -u)

  popd >/dev/null
done

echo "All peak calls & blacklist removal complete."






