# Pik3ca Single-Cell & Multi-Omics Analysis

This repository contains all the code and documentation used to analyze and visualize single-cell transcriptomics and 10X multi-omics datasets for our publication.

**Publication**

Currently, manuscript is on preparation.

> A DOI and link to the final publication will be provided here once available.

---

## Datasets

All raw sequencing datasets that support the findings of this study have been deposited in the NCBI Gene Expression Omnibus (GEO) under the following accession numbers.

- **10X Multiome** (GEO accession: [GSE282228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282228))  
- **10X scRNA-seq** (GEO accession: [GSE281982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281982))  

> **Note:** Processed Seurat objects and related data files are available upon request.

## Environment

All analyses were performed on **Ubuntu 18.04** using:

- **R** ≥ 4.4.2  
- **Conda** with dedicated environments for  
  - [pySCENIC](https://github.com/yurasong/Pik3ca_MG_analysis/blob/main/01_scRNA-seq/03_pySCENIC/conda_env.txt)  
  - [SCENIC+](https://github.com/yurasong/Pik3ca_MG_analysis/blob/main/03_SCENIC%2B/conda_environment_dependencies.txt)  

Package versions are listed in the **Methods** section of the manuscript and captured in the environment files.

## Installation

1. **Clone the repo**  
   ```bash
   git clone https://github.com/yourusername/Pik3ca_MG_analysis.git
   cd Pik3ca_MG_analysis

Each code file contains its own instructions for running the specific analyses.

---

## License

All rights are reserved by the authors.