#conda activate velocyto

velocyto run -b CTL/filtered_feature_bc_matrix/CTL_barcodes.tsv.gz -o CTL/ CTL/possorted_genome_bam.bam /home/audrey/reference_genomes/single-cell/refdata-gex-mm10-2020-A/genes/genes.gtf 

velocyto run -b ERPik/filtered_feature_bc_matrix/barcodes.tsv.gz -o ERPik/ ERPik/possorted_genome_bam.bam /home/audrey/reference_genomes/single-cell/refdata-gex-mm10-2020-A/genes/genes.gtf 

velocyto run -b KitPik/filtered_feature_bc_matrix/barcodes.tsv.gz -o KitPik/ KitPik/possorted_genome_bam.bam /home/audrey/reference_genomes/single-cell/refdata-gex-mm10-2020-A/genes/genes.gtf 
