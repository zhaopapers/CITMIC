# CITMIC
[![Current devel version: 0.3.0](https://img.shields.io/badge/devel%20version-0.3.0-blue.svg)](https://github.com/zhaopapers/CITMIC)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/09b138b2fa9242229f081cd180f6fc91)](https://app.codacy.com/gh/randef1ned/labyrinth/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](https://makeapullrequest.com)
[![Last commit](https://img.shields.io/github/last-commit/randef1ned/labyrinth.svg)](https://github.com/zhaopapers/CITMIC)
[![Code size](https://img.shields.io/github/languages/code-size/randef1ned/labyrinth.svg)](https://github.com/zhaopapers/CITMIC)

CITMIC is a method for estimate cell infiltration score

##Introduction

`CITMIC` can infer the cell infiltration of the TME by simultaneously measuring 86 different cell types, constructing an individualized cell-cell crosstalk network based on functional similarities between cells, and using only gene transcription data. This is a novel approach to estimate the relative cell infiltration levels, and which were shown to be superior to the current methods.

![A simple schema of the CITMIC](man/figures/info.jpg)

## A notice on operating system compatibility
- We recommend using as input the gene expression matrix normalized by log2(fpkm+1).
- If your computer has multiple CPU cores, we recommend that you select the number of cl.cores to apply to CITMIC based on `parallel::detectCores()`. He will greatly improve the efficiency of the function's operation.
- **R (≥ 4.2.0)**: We developed this R package using R version 4.3.x.
  

## Installation

Install `labyrinth` using:

``` r
install.packages(c('devtools', 'BiocManager'))
remotes::install_github("zhaopapers/CITMIC")
```

## Usage

Load the package using `library(CITMIC)`. We provide a vignette for the package that can be called using: `vignette("CITMIC")`. 



