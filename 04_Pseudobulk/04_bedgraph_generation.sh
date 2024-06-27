#!/bin/bash

##############################################################################
# This run has been proceeded on Ubuntu 18.04.
# All can be done on terminal. This is shell based.
# Scale factor is calculated as (normalized into 1M)/(rawcount) on each sample.
##############################################################################

bedtools genomecov -ibam HY_BC_ER+/K8Pik/paired_end_mapped.bam -bg -scale 11.71 > scaled_K8Pik_HY_BC_ER+.bedgraph
bedtools genomecov -ibam HY_BC_ER+/Klf5KO/paired_end_mapped.bam -bg -scale 33.40 > scaled_Klf5KO_HY_BC_ER+.bedgraph

bedtools genomecov -ibam HY_ER+_ER-/K8Pik/paired_end_mapped.bam -bg -scale 4.29 > scaled_K8Pik_HY_ER+_ER-.bedgraph
bedtools genomecov -ibam HY_ER+_ER-/Klf5KO/paired_end_mapped.bam -bg -scale 46.00 > scaled_Klf5KO_HY_ER+_ER-.bedgraph

bedtools genomecov -ibam Immature_BC/K8Pik/paired_end_mapped_K8Pik_Immature_BC.bam -bg -scale 1.49 > scaled_K8Pik_Immature_BC.bedgraph
bedtools genomecov -ibam Immature_BC/Klf5KO/paired_end_mapped_Klf5KO_ImmatureBC.bam -bg -scale 2.97 > scaled_Klf5KO_Immature_BC.bedgraph

bedtools genomecov -ibam LC_ER+/K8Pik/paired_end_mapped.bam -bg -scale 1.78 > scaled_K8Pik_LC_ER+.bedgraph
bedtools genomecov -ibam LC_ER+/Klf5KO/paired_end_mapped.bam -bg -scale 0.90 > scaled_Klf5KO_LC_ER+.bedgraph

bedtools genomecov -ibam LC_ER-/K8Pik/paired_end_mapped.bam -bg -scale 3.76 > scaled_K8Pik_LC_ER-.bedgraph
bedtools genomecov -ibam LC_ER-/Klf5KO/paired_end_mapped.bam -bg -scale 1.00 > scaled_Klf5KO_LC_ER-.bedgraph

bedtools genomecov -ibam K8PIk_only/Myoepith/paired_end_mapped.bam -bg -scale 4.25 > scaled_K8Pik_Myoepith.bedgraph
bedtools genomecov -ibam K8PIk_only/HY_BC_ER-/paired_end_mapped.bam -bg -scale 3.57 > scaled_K8Pik_HY_BC_ER-.bedgraph