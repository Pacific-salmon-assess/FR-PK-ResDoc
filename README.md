# Fraser-Pink-ResDoc
CSAS Research Document for Fraser Pink Salmon. This repository contains all the data, analysis, and csasdown code to fully reproduce the Research Document. We compile available data, fit a state-space spawner-recruit model, estimate biological benchmarks ($S_{MSY}$, $U_{MSY}$ and $S_{GEN}$), use a closed-loop simulation to project population dynamics forward under different harvest control rules (HCRs), then evaluate these dynamics and HCRs with how they perform in relation to conservation and fishery objectives.  
The root directory contains all files produced when running [csasdown](https://github.com/pbs-assess/csasdown)'s `draft("sr")` function. In addition to the csasdown files, 📁`analysis` contains:  
- 📁`data`  
  -📁`generated` data generated from analyses  
  -📁`raw ` data used in the analysis or to show the HCR  
    -📁`bio` biological data used to add information in the context section of the SR     
- 📁`R` `.R` files to fit the model (`fir-sr-stan.R`), run the forward simulation (`fwd-sim.R`), and create figures/tables (`inference-figs.R`)  
- 📁`Stan` Stan code to run the state-space spawner-recruit model  

To reproduce this document, fork or clone this repository, run the analysis, [fit the model](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/fit-sr-stan.R), run the [forward simulation](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/fwd-sim.R), [generate figures](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/inference-figs.R), and finally knit the [csasdown](https://github.com/pbs-assess/csasdown) generated [index](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/index.Rmd) document. Flag problems in [issues](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/issues) or reach out to dylan.glaser@dfo-mpo.gc.ca for help. 

This analysis depends on you loading: 
```{r env}
suppressWarnings({suppressMessages({
library(here)
library(tidyverse)
library(rosettafish)
library(csasdown)
library(kableExtra)
library(rstan)
library(gsl)
library(cowplot)
library(scales)
})})
pkgs <- sort(c("here", "tidyverse", "rosettafish", "csasdown", "kableExtra", "rstan", "gsl", "cowplot", "scales"))

devtools::session_info()$packages %>%
  dplyr::filter(package %in% pkgs) %>%
  dplyr::select(package, loadedversion, date) %>%
  dplyr::rename(Package = package, Version = loadedversion, Date = date) %>%
  csasdown::csas_table()
```
