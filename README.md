# Single-Cell and Multi-Omics Analysis of Pik3ca-Driven Mammary Epithelial Plasticity

This repository contains all the code and documentation used to analyze and visualize single-cell transcriptomics and 10X multi-omics datasets for our publication.

**Publication**

Currently, manuscript is on preparation.

> A DOI and link to the final publication will be provided here once available.

---

## Datasets

All raw sequencing datasets that support the findings of this study have been deposited in the NCBI Gene Expression Omnibus (GEO) under the following accession numbers.

- **10X Multiome** (GEO accession: [GSE282228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282228))  
- **10X scRNA-seq (CTL, ER Pik Yfp, Kit Pik Yfp)** (GEO accession: [GSE281982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281982))
- **10X scRNA-seq (ER Yfp-ER Pik Yfp, ER Pik Yfp - ER Pik Klf5)** (GEO accession: [GSE316159](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316159))
- **Human-relevance analyses**: [Wu_et_al_Cancer_epith_only.rds](https://www.dropbox.com/scl/fi/85ny19cwu73ugfvvi5qfs/Wu_et_al_Cancer_epith_only.rds?rlkey=2yw29ff6u53npgir3rpw4yndu&st=c4j4ff1e&dl=0)

> **Note:** Processed Seurat objects and related data files are available via [here](https://www.dropbox.com/scl/fo/2376enag74fwaewsz0uq5/AIHOMZih7uUJ1UGmmqeD2GM?rlkey=7oa4bxm2i4236iniko6ccstk9&st=ut0j445l&dl=0).

## Environment

All analyses were performed on **Ubuntu 20.04** and **MacOS 15.6 (24G84)** using:

- **R** ≥ 4.4.2  
- **Conda** with dedicated environments for  
  - [pySCENIC](https://github.com/yurasong/Pik3ca_MG_analysis/blob/main/01_scRNA-seq/04_pySCENIC/conda_env.txt)  
  - [SCENIC+](https://github.com/yurasong/Pik3ca_MG_analysis/blob/main/05_SCENIC%2B/conda_environment_dependencies.txt)  

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
| Cédric Blanpain     | **Corresponding author** (Study supervision) | [Cedric.Blanpain@ulb.be](mailto:Cedric.Blanpain@ulb.be) |
| Alexandra van Keymeulen | **Corresponding author** (Experimental) | [alexandra.van.keymeulen@ulb.be](mailto:alexandra.van.keymeulen@ulb.be) |
| Alejandro Sifrim    | **Corresponding author** (Bioinformatics) | [alejandro.sifrim@kuleuven.be](mailto:alejandro.sifrim@kuleuven.be) |
| Marco Fioramonti    | **First author** (Experimental)        | [ma.fioramonti@gmail.com](mailto:ma.fioramonti@gmail.com) |
| Yura Song           | **First author** (Bioinformatics)      | [yura.song@ulb.be](mailto:yura.song@ulb.be) |
| Brigida Novello     | **First author** (Experimental)        | [brigida.novello@ulb.be](mailto:brigida.novello@ulb.be) |
