# sFlt1-paper
Modeling and analysis code for sFLT1 trafficking manuscript (Gill et al., in progress)









## Matlab Code

### `01_ode_model/`

This folder contains code to simulate sFLT1 secretion as a system of 2 ordinary differential equations describing intracellular and extracellular sFLT1 pools.

To run simulations related to Figures 1 and S1, run the driver file `ode_driver.m`. All other scripts in the folder are called within this driver.

#### Driver run modes

- `"base_constitutive"`: run one constitutive secretion simulation with default parameter values (uses `sim_secr_ode.m`)
- `"base_on_off"`: run one simulation with equal periods of production and no production (uses `sim_on_off_ode.m`)
- `"base_pulse_chase"`: run one pulse-chase simulation with default parameter values (uses `sim_pulse_chase_ode.m`)
- `"extended_pulse_chase"`: run one pulse-chase simulation for longer time (uses `sim_pulse_chase_ode.m`)
- `"opt_1"`: run one nonlinear least squares optimization starting from default parameter values
- `"opt_100"`: run 100 optimizations from randomly sampled initial values

#### Other files

Each driver mode calls some combination of these files.

##### Setup scripts and helper functions

- `make_outdir.m`: Function to create a timestamped output directory
- `module_setup.m`: Adjust initial parameter values, system geometry, and time spans
- `module_hjk_opt_setup.m`: Adjust optimization options: cost function format, number of parameter sets, initial distributions of parameters, lower/upper bounds
- `module_init_hjk_arrays.m`: Initialize output containers for simulation data

##### Simulation functions

- `ode_eqns.m`: Equations used by `ode15s` solver to simulate sFLT1 secretion. Called by each of the `sim_` functions below.
- `sim_secr_ode.m`: Simulate constitutive sFLT1 production 
- `sim_on_off_ode.m`: Simulate a pulse of sFLT1 production without media change at end of pulse (just stop production).
- `sim_pulse_chase_ode.m`: Simulate a pulse of labeled sFLT1 production followed by chase (stop production and set extracellular sFLT1 to 0)

##### Helper functions, simulation modules

- `hjk_combo_cost.m`: Cost function to optimize trafficking parameters to simultaneously fit Jung et al., 2012 pulse-chase data, one Hornig et al., 2000 dataset, and control data from Kinghorn et al., 2023, within lsqnonlin using specified equations
- `module_constitutive_run_optimal.m`: Run constitutive simulation with optimal parameter values and report cost
- `module_jung_run_optimal.m`: Run pulse-chase simulation with optimal parameter values and report cost
