Recent Adaptation in Copadichromis mloto — Analysis Code Archive

This repository contains the reproducible code, workflows, and notebooks supporting the study on recent adaptation and fisheries-induced evolution in Copadichromis mloto populations from Lake Malombe and Lake Malawi.

Contents
```
Code_archive_mss/
├── Notebooks/ 
│   ├── recent_adaptation_Malombe_manuscript_code_pub.ipynb # Main manuscript analysis pipeline
│   ├── simulation_output_pub.ipynb # Explorations of coalescent simulation output
│   └── pyDeseq_pub.ipynb # DESEQ analyses
├── Snakemake/ # Workflow automations referenced in notebooks.
│   ├── Iterate_SNPeff # to get MAF corrected expectations of SNPeff catagories
│   ├── SFS_TD_FWH # Estimation of genome wide tajima's D and Fay and Wu's H from SFS
│   ├── Variantcalling # pipeline for variantcalling.
│   ├── fst # PBS calculations from VCF file
│   └── hscan # workflow for hscan calculation
├── Neutral_simulations/ # workflow automation and script to run demographic simulations to estimate neutral expectations for PBS and HSCAN.
├── TopGO/ # Contains Rscript for running topGO analysis
└── README.md
```
Usage

Publicly available datasets are not stored in this repository, as duplicating external resources is unnecessary and redundant. Versions or accession numbers are provided where relevant. VCF file for c.mloto snp data has been made available on Zenodo:10.5281/zenodo.17611130

Contact

For questions regarding the analysis or reproduction of results, please contact alex.hooftvanhuysduynen@uantwerpen.be


