# sFlt1-paper
Modeling and analysis code for sFLT1 trafficking manuscript (Gill et al., in progress)

## Prerequisites

Matlab: This code requires the Optimization Toolbox.

R: This code is written mainly in R markdown (Rmd) notebooks, which work best in the RStudio IDE. Required packages are listed in `helper_scripts/setup.R`.

## `01_ode_model/`

This folder contains code to simulate and analyze sFLT1 secretion as a system of 2 ordinary differential equations 
describing intracellular and extracellular sFLT1 pools. There is also code to optimize the model to experimental data.

### `01_ode_model/Matlab` 

#### Driver file: `ode_driver.m`

To run simulations related to Figures 1 and S1, run the driver file `ode_driver.m`. 
All other scripts in the folder are called within this driver.

##### Driver run modes

- `"base_constitutive"`: run one constitutive secretion simulation with default parameter values (uses `sim_secr_ode.m`)
- `"base_on_off"`: run one simulation with equal periods of production and no production (uses `sim_on_off_ode.m`)
- `"base_pulse_chase"`: run one pulse-chase simulation with default parameter values (uses `sim_pulse_chase_ode.m`)
- `"extended_pulse_chase"`: run one pulse-chase simulation for longer time (uses `sim_pulse_chase_ode.m`)
- `"opt_1"`: run one nonlinear least squares optimization starting from default parameter values
- `"opt_100"`: run 100 optimizations from randomly sampled initial values

##### Exporting simulation output data

To export output data for analysis in R, set `save_params = 1` and `save_data = 1`, 
which will save timestamped output files in `01_ode_model/Matlab/runs`. 
To run R scripts with these data, manually move desired output files to `saved-data/01_ode_model`.

#### Other files

Each driver mode calls some combination of these files.

##### Setup scripts and helper functions

- `make_outdir.m`: Function to create a timestamped output directory
- `module_hjk_opt_setup.m`: Adjust optimization options: cost function format, number of parameter sets, initial distributions of parameters, lower/upper bounds
- `module_init_hjk_arrays.m`: Initialize output containers for simulation data
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions and modules

- `hjk_combo_cost.m`: Cost function to optimize trafficking parameters to simultaneously fit Jung et al., 2012 pulse-chase data, one Hornig et al., 2000 dataset, and control data from Kinghorn et al., 2023, within lsqnonlin using specified equations
- `ode_eqns.m`: Equations used by `ode15s` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `sim_on_off_ode.m`: Simulate a pulse of sFLT1 production without media change at end of pulse (just stop production).
- `sim_pulse_chase_ode.m`: Simulate a pulse of labeled sFLT1 production followed by chase (stop production and set extracellular sFLT1 to 0)
- `sim_secr_ode.m`: Simulate constitutive sFLT1 production 

##### Additional modules

- `module_constitutive_run_optimal.m`: Run constitutive simulation with optimal parameter values and report cost
- `module_jung_run_optimal.m`: Run pulse-chase simulation with optimal parameter values and report cost

### `01_ode_model/R` 

IMPORTANT: To run R code, manually move desired simulation output files from timestamped `runs` output folders to `saved-data/01_ode_model`.

Any new datasets are saved to `saved-data/01_ode_model`. Figure outputs are saved to `saved-figs/01_ode_model`.

- `ode_toy_analysis.Rmd`: Uses output from the driver with mode "base_on_off" to produce Figure S1 panels.
- `ode_hjk_opt_import.Rmd`: Uses output from the driver with mode "opt_100" to import ODE optimization results.
- `ode_opt_workspace_2023-08-14.Rmd`: Uses imported ODE optimization results to produce Figure 1C-E panels.

## `02_dde_model`

This folder contains code to simulate sFLT1 secretion as a system of 2 delay differential equations 
describing intracellular and extracellular sFLT1 pools. 
There is a fixed time delay between production and either secretion or intracellular degradation. 
There is also code to optimize the model to experimental data.

### `02_dde_model/Matlab`

#### Driver file: `base_dde_driver.m`

To run simulations related to Figures (2, 4, 5, S2, S4, S5, S6, S7, S8), run the driver file `base_dde_driver.m`. All other scripts in the folder are called within this driver.

##### Driver run modes

- `"base_on_off"`: run one simulation with equal periods of production and no production
- `"vary_tau_on_off"`: run simulations as in `"base_on_off"` with a range of tau values
- `"base_model_opt"`: run optimization to experimental data 1000x (both pulse-chase and constitutive secretion)

##### Exporting simulation output

To export output data for analysis in R, set `save_params = 1` and `save_data = 1`, which will save timestamped output files in `02_dde_model/runs`. To run R scripts with these data, manually move output files to ________.

#### Other files

Each driver mode calls some combination of these files.

##### Setup scripts and helper functions

- `make_outdir.m`: Function to create a timestamped output directory
- `module_hjk_opt_setup.m`: Adjust optimization options: cost function format, number of parameter sets, initial distributions of parameters, lower/upper bounds
- `module_init_hjk_arrays.m`: Initialize output containers for simulation data
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions

- `dde_eqns.m`: Equations used by `dde23` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `hjk_combo_cost.m`: Cost function to optimize trafficking parameters to simultaneously fit Jung et al., 2012 pulse-chase data, one Hornig et al., 2000 dataset, and control data from Kinghorn et al., 2023, within lsqnonlin using specified equations
- `sim_on_off_dde.m`: Simulate a pulse of sFLT1 production without media change at end of pulse (just stop production).
- `sim_pulse_chase_dde.m`: Simulate a pulse of labeled sFLT1 production followed by chase (stop production and set extracellular sFLT1 to 0)
- `sim_secr_dde.m`: Simulate constitutive sFLT1 production 

##### Additional modules

- `module_constitutive_run_optimal.m`: Run constitutive simulation with optimal parameter values and report cost
- `module_jung_run_optimal.m`: Run pulse-chase simulation with optimal parameter values and report cost

### `02_dde_model/R`

IMPORTANT: To run R code, manually move desired simulation output files from timestamped `runs` output folders to `saved-data/02_dde_model`.

Any new datasets are saved to `saved-data/02_dde_model`. Figure outputs are saved to `saved-figs/02_dde_model`.

- `dde_toy_analysis.Rmd`: Uses output from the driver with mode "base_on_off" to produce Figure S2 panels (general solution of DDE model).
- `dde_vary_tau.Rmd`: Uses output from the driver with mode "vary_tau_on_off" to analyze the general solution of the DDE model at various tau levels (not in manuscript).
- `dde_opt_import.Rmd`: Uses output from the driver with mode "base_model_opt" to import DDE optimization results
- `dde_opt_preprocessing.Rmd`: Uses output from `dde_opt_import.Rmd`, preprocesses and filters, and produces Figure S4A,B 
- `dde_opt_workspace.Rmd`: Analyzes output from `dde_opt_preprocessing.Rmd` to produce these figures: 2B, 4, 5, S4C, S5, S6, S7, S8

## `03_candidate_models`

This folder contains code to optimize 8 different systems of ordinary and delay differential equations to experimental data. The equations include some combination of temporal decay in production rate, delayed intracellular clearance, and internalization of extracellular sFLT1.

### `03_candidate_models/Matlab`

#### Driver file: `driver_dde_candidate.m`

To run simulations related to Figures 3 and S3, run the driver file `driver_dde_candidate.m`. All other scripts in the folder are called within this driver.

##### Driver run modes

- `"hornig_basic_run"`: single run for troubleshooting 
- `"candidate_model_opt"`: multiple optimizations of each candidate model; adjust number of optimizations per model in file "module_hjk_opt_setup.m"

##### Exporting simulation outputs

To export output data for analysis in R, set `save_params = 1` and `save_data = 1`, which will save timestamped output files in `03_candidate_models/runs`. To run R scripts with these data, manually move output files to ________.

#### Other files

Each driver mode calls some combination of these files.

##### Setup scripts and helper functions

- `make_outdir.m`: Function to create a timestamped output directory
- `module_hjk_opt_setup.m`: Adjust optimization options: cost function format, number of parameter sets, initial distributions of parameters, lower/upper bounds
- `module_init_hjk_arrays.m`: Initialize output containers for simulation data
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions

- `dde_fun.m`: Equations used by `dde23` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `hjk_combo_cost.m`: Cost function to optimize trafficking parameters to simultaneously fit Jung et al., 2012 pulse-chase data, one Hornig et al., 2000 dataset, and control data from Kinghorn et al., 2023, within lsqnonlin using specified equations
- `sim_on_off_dde.m`: Simulate a pulse of sFLT1 production without media change at end of pulse (just stop production).
- `sim_pulse_chase_dde.m`: Simulate a pulse of labeled sFLT1 production followed by chase (stop production and set extracellular sFLT1 to 0)
- `sim_secr_dde.m`: Simulate constitutive sFLT1 production 

##### Additional modules

- `module_constitutive_run_optimal.m`: Run constitutive simulation with optimal parameter values and report cost
- `module_jung_run_optimal.m`: Run pulse-chase simulation with optimal parameter values and report cost

### `03_candidate_models/R`

- `candidate_models_import.Rmd`: Uses output from the driver with mode "candidate_model_opt" to import optimization results for 8 sets of DDEs
- `candidate_opt_preprocessing.Rmd`: Uses output from `candidate_models_import.Rmd`, preprocesses and filters
- `candidate_model_workspace.Rmd`: Analyzes output from `dde_opt_preprocessing.Rmd` to produce these figures: 3BC, S3

## `04_sensitivity`

This folder contains code that performs various types of sensitivity analysis on the optimized system of delay differential equations.

### `04_sensitivity/Matlab`

#### Driver file: `sens_analysis.m`

To run simulations related to Figures 5-8 and S9-S13, run the driver file `sens_analysis.m`. All other scripts in the folder are called within this driver.

##### Driver run modes

- `"base_const"`: constitutive secretion with base parameters
- `"local_sens_const"`: local (x%) sensitivity analysis of constitutive case parameters 
- `"sens_vary_inh_genchem"`: univariate sensitivity analysis - percent inhibition of individual parameters; change inhibition type with `genetic_chemical` setting
- `"sens_logscale_sqrt10"`: univariate sensitivity analysis - scale individual parameters by powers of sqrt(10)
- `"explore_c1"`: Vary alpha and beta while keeping c1 = alpha*beta constant
- `"explore_c2"`: Vary beta and gamma while keeping c2 = beta+gamma constant
- `"explore_c1c2"`: Vary alpha, beta, and gamma while keeping c1 = alpha*beta and c2 = beta+gamma constant
- `"sens_vary_base_inh_genchem"`: multivariate sensitivity analysis - from 5 different sets of initial {alpha, beta, gamma}, test percent inhibition of individual parameters

#### Other files

Each driver run mode calls some combination of these files

##### Setup scripts

- `make_outdir.m`: Function to create a timestamped output directory
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions

- `dde_eqns.m`: Equations used by `dde23` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `sim_secr_dde.m`: Simulate constitutive sFLT1 production 
- `sim_secr_dde_inh_t0.m`: Simulate constitutive sFLT1 production 

### `04_sensitivity/R`

All files here use output from "sFlt-model/04_sensitivity/sens_analysis.m" with the listed run modes. Files are listed in the order their figures appear in the manuscript.

- `local_sens_const.Rmd`: Uses output from run_mode = "local_sens_const" to produce Figure 6A 
- `powers_sqrt_10_constitutive.Rmd`: Uses output from run_mode = "sens_logscale_sqrt10" to produce Figures 6B and S9
- `explore_c1.Rmd`: Uses output from run_mode = "explore_c1" to produce Figure S10AB
- `explore_c2.Rmd`: Uses output from run_mode = "explore_c2" to produce Figure S10CD
- `frac_inh_chemical.Rmd`: Uses output from run_mode = "sens_vary_inh_genchem" with genetic_chemical = "chemical" to produce Figure 7BCD, S11A
- `frac_inh_genetic.Rmd`: Uses output from run_mode = "sens_vary_inh_genchem" with genetic_chemical = "genetic" to produce Figure 7EFG, S11B
- `diffbase_frac_inh_chemical.Rmd`: Uses output from run_mode = "sens_vary_base_inh_genchem" with genetic_chemical = "chemical" to produce Figure 8A, S12, S13A
- `diffbase_frac_inh_genetic.Rmd`: Uses output from run_mode = "sens_vary_base_inh_genchem" with genetic_chemical = "genetic" to produce Figure 8B, S12 S13B

## saved-data

### Initial content on Github

Initially contains experimental datasets from listed references and metadata for some optimizations.

#### Experimental datasets

- `hornig.mat`: Matlab struct containing experimental data from Hornig et al., 2000
- `hornig.rda`: R data frame containing experimental data from Hornig et al., 2000
- `jung_1d.rda`: Matlab struct containing experimental data from Jung et al., 2012
- `jung_v2.mat`: R data frame containing experimental data from Jung et al., 2012
- `kinghorn_ctrl_summary.mat`: Matlab struct containing experimental control data from Kinghorn et al., 2024
- `inh_data/inh_data_wrangled.csv`: Table with experimental data for chemical inhibitors from Kinghorn et al., 2024
- `inh_data/kd_wrangled.csv`: Table with experimental data for genetic knockdowns from Kinghorn et al., 2024

#### Metadata files

- `02_dde_model/hjk_dde_opt_meta.csv`: metadata table used by `02_dde_model/R/dde_opt_import.Rmd` to correctly format imported simulation output data
- `03_candidate_models/hjk_opt_candidates_meta.csv`: metadata table used by `03_candidate_models/R/candidate_models_import.Rmd` to correctly format imported simulation output data

#### Simulation output file generation

Simulation output files are not available through GitHub and must be generated using the Matlab scripts above.

R scripts will look in this directory for simulation output data. Manually move desired output data from the timestamped `runs` directory to this directory. Some output files may require renaming to match file names in R scripts.

## References

### Manuscript

bioRxiv: TBA

### Experimental data

Hornig, C. et al. Release and Complex Formation of Soluble VEGFR-1 from Endothelial Cells and Biological Fluids. Lab. Invest. 80, 443–454 (2000).

Jung, J.-J. et al. Secretion of Soluble Vascular Endothelial Growth Factor Receptor 1 (sVEGFR1/sFlt1) Requires Arf1, Arf6, and Rab11 GTPases. PLoS ONE 7, e44572 (2012).

Kinghorn, K. et al. A defined clathrin-mediated trafficking pathway regulates sFLT1/VEGFR1 secretion from endothelial cells. Angiogenesis 27, 67–89 (2024).

