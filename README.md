# ALIES
# ALIES: Arm-Level Immune Exclusion Score

## Overview

This repository contains the complete reproducible analysis pipeline for
the ALIES pan-cancer study, covering arm-level CNA correlation analysis,
LASSO score construction, survival analysis, and IMvigor210 external
validation.

## Requirements

- R >= 4.3.0
- Bioconductor >= 3.17

Install R packages by running the top of `R/ALIES_analysis.R` — all
packages are installed automatically on first run.

## Input Data

Download the following files and place them in your working directory:

| File | Source | How to download |
|------|--------|----------------|
| `PANCAN_Aneuploidy.xlsx` | Taylor et al. 2018 | Supplementary Table 2 of doi:10.1016/j.ccell.2018.03.007 |
| `TCGA_clinical_data.tsv.gz` | GDC Data Portal | https://portal.gdc.cancer.gov/ → Repository → Clinical |
| `TCGA_expression_matrix.tsv.gz` | GDC Data Portal | Pan-Cancer RNA-seq FPKM-UQ |
| `TCGA_mutation.maf.gz` | GDC Data Portal | MC3 MAF file |
| IMvigor210CoreBiologies | Genentech | http://research-pub.gene.com/IMvigor210CoreBiologies/ |

## How to Run

```r
# 1. Set your working directory to where the data files are
setwd("/path/to/your/data/folder")

# 2. Source the analysis script
source("R/ALIES_analysis.R")
```

All outputs are written to `outputs/tables/` and `outputs/figures/`.
Expected runtime: approximately 2–4 hours on a standard laptop
(immune deconvolution is the slowest step).

## Repository Contents
