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

For our human-relevance analyses (Wu et al., Nat Genet 2021), we downloaded the processed Seurat object used in the paper. You can retrieve it here: [Wu_et_al_Cancer_epith_only.rds](https://www.dropbox.com/scl/fi/85ny19cwu73ugfvvi5qfs/Wu_et_al_Cancer_epith_only.rds?rlkey=2yw29ff6u53npgir3rpw4yndu&st=c4j4ff1e&dl=0)

> **Note:** Processed Seurat objects and related data files are available upon request.

## Environment

All analyses were performed on **Ubuntu 20.04** and **MacOS 15.6 (24G84)** using:

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

---

## Contact

For requests of processed Seurat objects, related data files and protocols, please contact:

| Name                | Role                                   | Email |
|---------------------|----------------------------------------|-------|
| Cédric Blanpain     | **Corresponding author**               | [Cedric.Blanpain@ulb.be](mailto:Cedric.Blanpain@ulb.be) |
| Alexandra van Keymeulen | **Corresponding author** (Experimental) | [alexandra.van.keymeulen@ulb.be](mailto:alexandra.van.keymeulen@ulb.be) |
| Alejandro Sifrim    | **Corresponding author** (Bioinformatics) | [alejandro.sifrim@kuleuven.be](mailto:alejandro.sifrim@kuleuven.be) |
| Marco Fioramonti    | **First author** (Experimental)        | [ma.fioramonti@gmail.com](mailto:ma.fioramonti@gmail.com) |
| Yura Song           | **First author** (Bioinformatics)      | [yura.song@ulb.be](mailto:yura.song@ulb.be) |
| Brigida Novello     | **First author** (Experimental)        | [brigida.novello@ulb.be](mailto:brigida.novello@ulb.be) |