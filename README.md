# sFlt1-paper
Modeling and analysis code for sFLT1 trafficking manuscript (Gill et al., in progress)

## Prerequisites

Matlab: This code requires the Optimization Toolbox.

R: This code is written in R markdown (Rmd) notebooks, which work best in the RStudio IDE. Required packages are listed in `helper_scripts/setup.R`.

## `01_ode_model/`

This folder contains code to simulate and analyze sFLT1 secretion as a system of 2 ordinary differential equations 
describing intracellular and extracellular sFLT1 pools. There is also code to optimize the model to experimental data.

### Matlab 

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

#### R 

Any new datasets are saved to `saved-data/01_ode_model`. Figure outputs are saved to `saved-figs/01_ode_model`.

- `ode_toy_analysis.Rmd`: Uses output from the driver with mode "base_on_off" to produce Figure S1 panels.
- `ode_hjk_opt_import.Rmd`: Uses output from the driver with mode "opt_100" to import ODE optimization results.
- `ode_opt_workspace_2023-08-14.Rmd`: Uses imported ODE optimization results to produce Figure 1C-E panels.

### `02_dde_model

This folder contains code to simulate sFLT1 secretion as a system of 2 delay differential equations 
describing intracellular and extracellular sFLT1 pools. 
There is a fixed time delay between production and either secretion or intracellular degradation. 
There is also code to optimize the model to experimental data.

#### Driver file: `base_dde_driver.m`

To run simulations related to Figures ___ and ___, run the driver file `base_dde_driver.m`. All other scripts in the folder are called within this driver.

##### Driver run modes

- `"base_on_off"`: run one simulation with equal periods of production and no production
- `"vary_tau_on_off"`: run simulations as in `"base_on_off"` with a range of tau values
- `"base_model_opt"`: run optimization to experimental data 1000x (both pulse-chase and constitutive secretion)

##### Exporting simulation output

To export output data for analysis in R, set `save_params = 1` and `save_data = 1`, which will save timestamped output files in `01_ode_model/runs`. To run R scripts with these data, manually move output files to ________.

#### Other files

Each driver mode calls some combination of these files.

##### Setup scripts and helper functions

These files are analogous to the similarly named files in `01_ode_model`.

- `make_outdir.m`: Function to create a timestamped output directory
- `module_hjk_opt_setup.m`: Adjust optimization options: cost function format, number of parameter sets, initial distributions of parameters, lower/upper bounds
- `module_init_hjk_arrays.m`: Initialize output containers for simulation data
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions

These files are analogous to the similarly named files in `01_ode_model`.

- `dde_eqns.m`: Equations used by `dde23` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `hjk_combo_cost.m`: Cost function to optimize trafficking parameters to simultaneously fit Jung et al., 2012 pulse-chase data, one Hornig et al., 2000 dataset, and control data from Kinghorn et al., 2023, within lsqnonlin using specified equations
- `sim_on_off_dde.m`: Simulate a pulse of sFLT1 production without media change at end of pulse (just stop production).
- `sim_pulse_chase_dde.m`: Simulate a pulse of labeled sFLT1 production followed by chase (stop production and set extracellular sFLT1 to 0)
- `sim_secr_dde.m`: Simulate constitutive sFLT1 production 

##### Additional modules

These files are analogous to the similarly named files in `01_ode_model`.

- `module_constitutive_run_optimal.m`: Run constitutive simulation with optimal parameter values and report cost
- `module_jung_run_optimal.m`: Run pulse-chase simulation with optimal parameter values and report cost

### `03_candidate_models`

This folder contains code to optimize 8 different systems of ordinary and delay differential equations to experimental data. The equations include some combination of temporal decay in production rate, delayed intracellular clearance, and internalization of extracellular sFLT1.

#### Driver file: `driver_dde_candidate.m`

To run simulations related to Figures ___ and ___, run the driver file `driver_dde_candidate.m`. All other scripts in the folder are called within this driver.

##### Driver run modes

- `"hornig_basic_run"`: single run for troubleshooting 
- `"candidate_model_opt"`: multiple optimizations of each candidate model; adjust number of optimizations per model in file "module_hjk_opt_setup.m"

##### Exporting simulation outputs

To export output data for analysis in R, set `save_params = 1` and `save_data = 1`, which will save timestamped output files in `01_ode_model/runs`. To run R scripts with these data, manually move output files to ________.

#### Other files

Each driver mode calls some combination of these files.

##### Setup scripts and helper functions

These files are analogous to the similarly named files in previous modules.

- `make_outdir.m`: Function to create a timestamped output directory
- `module_hjk_opt_setup.m`: Adjust optimization options: cost function format, number of parameter sets, initial distributions of parameters, lower/upper bounds
- `module_init_hjk_arrays.m`: Initialize output containers for simulation data
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions

These files are analogous to the similarly named files in previous modules.

- `dde_fun.m`: Equations used by `dde23` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `hjk_combo_cost.m`: Cost function to optimize trafficking parameters to simultaneously fit Jung et al., 2012 pulse-chase data, one Hornig et al., 2000 dataset, and control data from Kinghorn et al., 2023, within lsqnonlin using specified equations
- `sim_on_off_dde.m`: Simulate a pulse of sFLT1 production without media change at end of pulse (just stop production).
- `sim_pulse_chase_dde.m`: Simulate a pulse of labeled sFLT1 production followed by chase (stop production and set extracellular sFLT1 to 0)
- `sim_secr_dde.m`: Simulate constitutive sFLT1 production 

##### Additional modules

These files are analogous to the similarly named files in previous modules.

- `module_constitutive_run_optimal.m`: Run constitutive simulation with optimal parameter values and report cost
- `module_jung_run_optimal.m`: Run pulse-chase simulation with optimal parameter values and report cost

### `04_sensitivity`

This folder contains code that performs various types of sensitivity analysis on the optimized system of delay differential equations.

#### Driver file: `sens_analysis.m`

To run simulations related to Figures ___ and ___, run the driver file `sens_analysis.m`. All other scripts in the folder are called within this driver.

##### Driver run modes

- `"base_const"`: constitutive secretion with base parameters
- `"local_sens_const"`: local (x%) sensitivity analysis of constitutive case parameters 
- `"sens_vary_inh_genchem"`: univariate sensitivity analysis - percent inhibition of individual parameters
- `"sens_logscale_sqrt10"`: univariate sensitivity analysis - scale individual parameters by powers of sqrt(10)
- `"explore_c1"`: Vary alpha and beta while keeping c1 = alpha*beta constant
- `"explore_c2"`: Vary beta and gamma while keeping c2 = beta+gamma constant
- `"explore_c1c2"`: Vary alpha, beta, and gamma while keeping c1 = alpha*beta and c2 = beta+gamma constant
- `"sens_vary_base_inh_genchem"`: multivariate sensitivity analysis - from 5 different sets of initial {alpha, beta, gamma}, test percent inhibition of individual parameters

#### Other files

Each driver run mode calls some combination of these files

##### Setup scripts

These files are analogous to the similarly named files in previous modules.

- `make_outdir.m`: Function to create a timestamped output directory
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans

##### Simulation functions

- `dde_eqns.m`: Equations used by `dde23` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `sim_secr_dde.m`: Simulate constitutive sFLT1 production 
- `sim_secr_dde_inh_t0.m`: Simulate constitutive sFLT1 production 

## R code

IMPORTANT: To run R code, manually move desired simulation output files from timestamped `runs` output folders to ...

