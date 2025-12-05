# brevi_etal_prostateCancer  
**Reanalysis of public prostate cancer transcriptomic datasets**  

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17832151.svg)](https://doi.org/10.5281/zenodo.17832151)

This repository provides scripts to re-analyse publicly available array, bulk RNA-seq, and single-cell RNA-seq datasets from prostate cancer studies.  

---

## 📄 Table of Contents

- [Datasets & Scripts](#datasets--scripts)  
- [Requirements](#requirements)  
- [Usage](#usage)  
- [License](#license)  
- [Contact](#contact)  
- [Acknowledgements](#acknowledgements)  

---

## Datasets & Scripts

### Single-cell RNA-seq  
| Study | Script | Data source / Notes |
|---|---|---|
| Cheng et al., *Eur Urol*, 2022 | `scRNA_analysis_Chengetal.R` | Data obtained directly from the authors — see https://doi.org/10.1016/j.eururo.2021.12.039 |
| Song et al., *Nat Commun*, 2022 | `scRNA_analysis_Songetal.R` | GEO: [GSE176031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176031) — see https://pmc.ncbi.nlm.nih.gov/articles/PMC8748675/ |
| Zaidi et al., *PNAS*, 2024 | `scRNA_analysis_Zaidietal.R` | GEO: [GSE264573](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264573) — see https://www.pnas.org/doi/10.1073/pnas.2322203121 |

### Bulk RNA-seq & Microarray  
| Study | Script | Data source |
|---|---|---|
| Beltran et al., *Nat Med*, 2016 | `RNAseq_analysis_Beltranetal.R` | dbGaP: phs000909.v.p1 — see https://www.nature.com/articles/nm.4045 |
| Jachetti et al., *Cancer Res*, 2015 | `array_analysis_Jachettietal.R` | GEO: [GSE65502](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65502) — see https://aacrjournals.org/cancerres/article/75/10/2095/599799 |
| Abida et al., *PNAS*, 2019 | `RNAseq_analysis_Abida.R` | Data kindly provided by A.Chinnaiyan (arul@umich.edu) and M.Cieslik (mcieslik@med.umich.edu) upon request — see https://doi.org/10.1073/pnas.1902651116 |

---

## Requirements  

- R (≥ 4.2)  
- Recommended R packages (depending on the script):  
  - `Seurat`, `SingleCellExperiment`, `DESeq2`, `limma`, `tidyverse`,  `dlpry`,  `openxlsx`, `ggplot2`


---

## Usage  

Each script can be run independently. Examples:

```bash
Rscript scRNA_analysis_Zaidietal.R  
Rscript RNAseq_analysis_Beltranetal.R  
```

To use you need to update the path of input data once retrived. 


Input file with tested signatures are provided in the folder `input/`.


## Contributors

- Anna Sofia Tascini (tascini.annasofia@unisr.it)

- Aurora Maurizio (maurizio.aurora@hsr.it)
