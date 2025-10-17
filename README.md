# Estimating Precautionary Approach reference points and assessing consequences of harvest control rules for Fraser River Pink salmon 

[CSAS Research Document for Fraser Pink Salmon](https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-5/Fs70-5-2024-063-eng.pdf). This repository contains all the data, analysis, and csasdown code to fully reproduce the Research Document. We compile available data, fit a state-space spawner-recruit model, estimate biological benchmarks ($S_{MSY}$, $U_{MSY}$ and $S_{GEN}$), use a closed-loop forward simulation conditioned on recent population dynamics to project the stock forward under different harvest control rules (HCRs). Then we evaluate HCRs by how they perform in relation to conservation and fishery objectives.  

The root directory contains all files produced when running [csasdown](https://github.com/pbs-assess/csasdown)'s `draft("sr")` function. In addition to the csasdown files, 📁`analysis` contains:  
- 📁`data`  
  -📁`generated` data generated from analyses  
  -📁`raw ` data used in the analysis or to show the HCR with subfolders:  📁`bio` containing biological data used to add information in the context section and 📁`hatchery` with hatchery supplementation data from SEP       
- 📁`R` `.R` files to fit the model (`fir-sr-stan.R`), run the forward simulation (`fwd-sim.R`), create figures/tables (`inference-figs.R`) and some helper functions  
- 📁`Stan` Stan code to run the state-space spawner-recruit models
- [The model validation supplement](https://pacific-salmon-assess.github.io/FR-PK-ResDoc/Supplement-model-check.html)

To reproduce this document: 
- fork or clone this repository
- make sure you have all [packages](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/10-appendix.Rmd#L25) and software (i.e. Stan) installed
- run the analysis (in *this* order) by:
  -  [fitting the model](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/fit-sr-stan.R),
  -  [running the forward simulation](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/fwd-sim.R),
  -  [generating the figures](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/analysis/R/inference-figs.R),
  -  and finally knit the [index.Rmd](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/blob/main/index.Rmd) file to generate the document.
- Flag problems relating to the analysis and structure of this repo in [issues](https://github.com/Pacific-salmon-assess/FR-PK-ResDoc/issues) or reach out to dylan.glaser@dfo-mpo.gc.ca for help. 

Citation: 

Glaser, D.M., Connors, B.M., Dionne, K., and Huang, A.M. 2025. Estimating Precautionary 
  Approach Reference Points and Assessing Consequences of Harvest Control Rules for 
  Fraser River Pink Salmon (*Oncorhynchus gorbuscha*). DFO Can. Sci. Advis. Sec. Res. Doc. 
  2024/063. iv + 31 p.