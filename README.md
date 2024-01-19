# Fraser-Pink-SR
CSAS Science Response (SR) for Fraser Pink Salmon. This repository contains all the data, analysis, and csasdown code to fully reproduce the Science Response. We compile available data, fit a state-space spawner-recruit model, estimate biological benchmarks ($S_{MSY}$, $U_{MSY}$ and $S_{GEN}$), use a closed-loop simulation to project population dynamics forward under different harvest control rules (HCRs), then evaluate these dynamics and HCRs with how they relate to estimated biological benchmarks and operational control points (OCPs).  
The root directory contains all files produced when running [csasdown](https://github.com/pbs-assess/csasdown)'s `draft("sr")` function. In addition to the csasdown files, 📁`analysis` contains:  
- 📁`data`  
  -📁`generated` data generated from analyses  
  -📁`raw ` data used in the analysis or to show the HCR  
    -📁`bio` biological data used to add information in the context section of the SR     
- 📁`R` `.R` files to fit the model (`fir-sr-stan.R`), run the forward simulation (`fwd-sim.R`), and create figures/tables (`inference-figs.R`)  
- 📁`Stan` Stan code to run the state-space spawner-recruit model  

To re-run the analysis, fit the model, run the forward simulation, then look at the figures and tables.  

**depends on** -add dependencies
