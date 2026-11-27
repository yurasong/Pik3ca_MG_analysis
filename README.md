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
- **Human-relevance analyses**: Publicly available data from [Wu_et_al., *Nat Genet*, 2021](https://www.nature.com/articles/s41588-021-00911-1)

## Environment

We performed the analyses on Ubuntu 20.04 and macOS 15.6 (24G84), using R ≥ 4.4.2 for the R analyses. We used separate Conda environments for [pySCENIC](./01_scRNA-seq/01-01_CTL-ERPik-KitPik/06_pySCENIC/environment_SCENIC.yml) and [SCENIC+](./02_multiome/03_SCENIC+/environment_SCENICplus.yml). Both YAML files specify package builds for 64-bit Linux (x86_64).

Additional software versions are listed in the Methods section of the manuscript.

## Installation

Clone the repository and enter its root directory:

```bash
git clone https://github.com/yurasong/Pik3ca_MG_analysis.git
cd Pik3ca_MG_analysis
```

On 64-bit Linux with Conda installed, follow the setup for the analysis you want to run.

### pySCENIC

```bash
conda env create -f 01_scRNA-seq/01-01_CTL-ERPik-KitPik/06_pySCENIC/environment_SCENIC.yml
conda activate pyscenic
pyscenic --help
```

### SCENIC+

The installation uses the dependency versions and build tools specified in the YAML file. The two `PIP_` settings apply only to the environment creation command.

```bash
PIP_NO_DEPS=1 PIP_NO_BUILD_ISOLATION=0 conda env create -f 02_multiome/03_SCENIC+/environment_SCENICplus.yml
conda activate scenicplus
python -c "from scenicplus.scenicplus_class import SCENICPLUS; print('SCENIC+ import successful')"
```

Each analysis script contains instructions for running the corresponding analysis.

---

## License

The analysis code in this repository is available under the [MIT License](LICENSE).

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
