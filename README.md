# Gametocytogenesis

<img src="https://github.com/andyrussell/Gametocytogenesis/blob/master/GCSKO_logo.jpg?raw=true" width="500">

In this project, we explore the process of gametocytogenesis in *Plasmodium berghei*.

Please contact ar19@sanger.ac.uk with any questions, or if you would like to collaborate.


## Reproducing our analysis

This work was conducted using R 4.05.

You can install necessary packages by installing `renv` and running `renv::restore()` with this directory as your working directory.

### Single-cell data
- 1. Run `quality_control` notebooks to generate QCed data files
    - `GCSKO_SS2_QC.Rmd`
    - `GCSKO_10X_QC.Rmd`
    - `GCSKO_genotyping_coverage_plots.Rmd`*
- 2. Run `analysis` notebooks
    - A. `merge`
    - B. `coexpression`
    - C. `sex_branch_analysis`
    - D. `sex_ratio`
    - E. `misc_code` 

*optional
