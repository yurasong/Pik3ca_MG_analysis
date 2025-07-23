#!/bin/bash

##########################################################################################
# File: 01_proceed_bamfile
#
# Description:
#   Filters ATAC-seq BAM files by cell barcodes and cleans contigs/reads per sample:
#     1) Uses `sinto filterbarcodes` to extract per-cell BAMs based on metadata
#     2) Strips non-standard contigs (keeps chr1–22, X, Y)
#     3) Merges, sorts, and indexes cleaned BAMs
#     4) Filters for mapped and paired reads
#
# Inputs:
#   - ATAC BAMs:
#       K8Pik_atac_possorted_bam.bam
#       K8Pik_Klf5KO_atac_possorted_bam.bam
#       CTL_atac_possorted_bam.bam
#   - Metadata files (tab-delimited):
#       00_input/metadata_K8Pik_noFibro_multiome.txt
#       00_input/metadata_K8Pik_Klf5KO_noFibro_multiome.txt
#       00_input/metadata_CTL_multiome.txt
#
# Outputs:
#   - `<sample>/<cell>.bam` (filtered by barcode)
#   - `<prefix>_noContig_sorted.bam` + index
#   - `<prefix>_mapped.bam`
#   - `<prefix>_paired.bam`
#
# Dependencies:
#   - bash
#   - samtools
#   - sinto (https://github.com/mojaveazure/sinto)
#   - coreutils (mktemp, mkdir, basename)
#
##########################################################################################

#----------------------------------------
# CONFIGURATION
#----------------------------------------
bams=( 
  "K8Pik_atac_possorted_bam.bam"
  "K8Pik_Klf5KO_atac_possorted_bam.bam"
  "CTL_atac_possorted_bam.bam"
)

metas=( 
  "00_input/metadata_K8Pik_noFibro_multiome.txt"
  "00_input/metadata_K8Pik_Klf5KO_noFibro_multiome.txt"
  "00_input/metadata_CTL_multiome.txt"
)

contigs=(chr{1..22} chrX chrY)

#----------------------------------------
# FUNCTIONS
#----------------------------------------
filter_barcodes() {
  local bam_file=$1 meta_file=$2
  local sample="${bam_file%%_atac_possorted_bam.bam}"
  echo "▶ Filtering ${bam_file} → ${sample}/"
  mkdir -p "$sample"
  sinto filterbarcodes \
    -b "$bam_file" \
    -c "$meta_file" \
    --outdir "$sample"
}

process_bam() {
  local in_bam=$1
  local prefix="${in_bam%.bam}"
  local tmpdir; tmpdir=$(mktemp -d)

  echo "▶ Stripping contigs from ${in_bam}"
  for c in "${contigs[@]}"; do
    samtools view -bh "$in_bam" "$c" > "${tmpdir}/${c}.bam"
  done

  echo "▶ Merging → ${prefix}_noContig.bam"
  samtools merge -O BAM "${prefix}_noContig.bam" "${tmpdir}"/*.bam
  rm -rf "$tmpdir"

  echo "▶ Sorting → ${prefix}_noContig_sorted.bam"
  samtools sort -O BAM -T "${prefix}_tmp" -o "${prefix}_noContig_sorted.bam" "${prefix}_noContig.bam"
  samtools index "${prefix}_noContig_sorted.bam"

  echo "▶ Filtering reads → mapped & paired"
  samtools view -b -F 4 "${prefix}_noContig_sorted.bam" > "${prefix}_mapped.bam"
  samtools view -b -f 1 "${prefix}_mapped.bam"      > "${prefix}_paired.bam"

  echo "✔ Done with ${in_bam}"
}

#----------------------------------------
# MAIN
#----------------------------------------

# 1) per-sample barcode filtering
for i in "${!bams[@]}"; do
  filter_barcodes "${bams[i]}" "${metas[i]}"
done
echo "✅ Barcode filtering complete."

# 2) per-cell-type contig & read filtering
for meta in "${metas[@]}"; do
  sample=$(basename "$meta")
  sample=${sample#metadata_}
  sample=${sample%_multiome.txt}

  echo "=== Processing sample: ${sample} ==="
  pushd "$sample" >/dev/null || { echo "No dir $sample, skipping"; continue; }

  # grab unique cell-type names from 2nd column
  while read -r cell; do
    bam="${cell}.bam"
    if [[ -f "$bam" ]]; then
      process_bam "$bam"
    else
      echo "⚠️  Missing ${bam}, skipping"
    fi
  done < <(cut -f2 "../$meta" | sort -u)

  popd >/dev/null
done

echo "All contig-stripping & read-filtering done."