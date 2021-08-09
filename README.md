# Gametocytogenesis

<img src="https://github.com/andyrussell/Gametocytogenesis/blob/master/GCSKO_logo.jpg?raw=true" width="500">

In this project, we explore the process of gametocytogenesis in *Plasmodium berghei*.

Source code to reproduce the analysis of the paper:

>**Regulators of male and female sexual development critical for transmission of a malaria parasite**
<br>Andrew J C Russell, Theo Sanderson, Ellen Bushell, Arthur M Talman, Burcu Anar, Gareth Girling, Mirjam Hunziker, Robyn S Kent, Tom Metcalf, Ruddy Montandon, Vikash Pandey, A Brett Roberts, Claire Sayers, Frank Schwach, Julian C Rayner, Thierry Voet, Katarzyna K Modrzynska, Andrew P. Waters, Mara K N Lawniczak, Oliver Billker<br>
DOI: https://doi.org/10.1101/2021.08.04.455056 

Please contact ar19@sanger.ac.uk with any questions, or if you would like to collaborate.


## Basic steps

### Single-cell data
- 1. Run `quality_control` notebooks to generate QCed data files
    - `GCSKO_SS2_QC.Rmd`
    - `GCSKO_10X_QC.Rmd`
    - `GCSKO_genotyping_coverage_plots.Rmd`*
- 2. Run `analysis` notebooks
    - A. `merge_wt`
    - B. `sex_branch_analysis_wt`
    - C. `branch_analysis`
    - D. `merge`
    - E. `sex_branch_analysis`
    - F. `sex_ratio`

*optional - slow run time

### Screen data

located in `barseq` directory.

### Bulk RNA-seq data

located in `bulk_rna-seq` directory.
