# CITMIC
[![Current devel version: 0.1.2](https://img.shields.io/badge/devel%20version-0.1.2-blue.svg)](https://github.com/zhaopapers/CITMIC)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2ba2ad32650d469588a16de5ae2a5ed1)](https://app.codacy.com/gh/zhaopapers/CITMIC/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

CITMIC is a method for estimate cell infiltration score

## Introduction

`CITMIC` can infer the cell infiltration of the TME by simultaneously measuring 86 different cell types, constructing an individualized cell-cell crosstalk network based on functional similarities between cells, and using only gene transcription data. This is a novel approach to estimate the relative cell infiltration levels, and which were shown to be superior to the current methods.

![A simple schema of the CITMIC](man/figures/info.jpg)

## A notice on operating system compatibility
- We recommend using as input the gene expression matrix normalized by log2(fpkm+1).
- If your computer has multiple CPU cores, we recommend that you select the number of cl.cores to apply to CITMIC based on `parallel::detectCores()`. He will greatly improve the efficiency of the function's operation.
- **R (≥ 4.2.0)**: We developed this R package using R version 4.3.x.
  

## Installation

Install `CITMIC` using:

``` r
install.packages(c('devtools', 'BiocManager'))
remotes::install_github("zhaopapers/CITMIC")
```

## Usage

Load the package using `library(CITMIC)`. We provide a vignette for the package that can be called using: `vignette("CITMIC")`. 
Alternatively, you can view the online version on [CRAN](doc/labyrinth_knit.md), The examples I provided would take several minutes to run on a normal desktop computer. Basically that is all you have to know.

## Citation
Those codes and the CITMIC package are intended for research use only. 

If you use CITMIC or these codes in your publication, please cite the paper: 

X. Zhao, J. Wu, J. Lai, B. Pan, M. Ji, X. Li, Y. He, J. Han, CITMIC: Comprehensive Estimation of Cell Infiltration in Tumor Microenvironment based on Individualized Intercellular Crosstalk. Adv. Sci. 2025, 12, 2408007. https://doi.org/10.1002/advs.202408007
