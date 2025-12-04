# Reanalysis of Public Prostate Cancer Transcriptomic Datasets

This repository contains scripts used to reanalyse public array, bulk RNA-seq, and single-cell RNA-seq datasets from prostate cancer studies.

---

## 📘 Single-cell RNA-seq

We reanalysed three single-cell RNA-seq datasets:

### • Cheng et al., *Eur Urol*, 2022  
Script: `scRNA_analysis_Chengetal.R`  
Input data obtained directly from the authors.

### • Song et al., *Nat Commun*, 2022  
Script: `scRNA_analysis_Songetal.R`  
Data retrieved from GEO: [GSE176031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176031)  
Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8748675/

### • Zaidi et al., *PNAS*, 2024  
Script: `scRNA_analysis_Zaidietal.R`  
Data retrieved from GEO: [GSE264573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264573)  
Paper: https://www.pnas.org/doi/10.1073/pnas.2322203121

---

## 📗 Bulk RNA-seq and Microarray

### • Beltran et al., *Nat Med*, 2016  
Script: `RNAseq_analysis_Beltranetal.R`  
Data obtained from dbGaP: **phs000909.v.p1**

### • Jachetti et al., *Cancer Res*, 2015  
Script: `array_analysis_Jachettietal.R`  
Data retrieved from GEO: [GSE65502](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65502)  
Paper: https://aacrjournals.org/cancerres/article/75/10/2095/599799

---

## 📦 Requirements

- R ≥ 4.2  
- Recommended packages: `Seurat`, `SingleCellExperiment`, `DESeq2`, `limma`, `tidyverse`, etc.  
(You can adjust this section based on your actual scripts.)

---

## ▶️ Usage

Each script can be run independently.  
Example:

```bash
Rscript scRNA_analysis_Zaidietal.R
