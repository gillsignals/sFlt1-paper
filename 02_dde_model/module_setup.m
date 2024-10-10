% set rng seed
rng(072823)

% if anything will be saved, create timestamped output directory
if (save_params || save_data || save_figs)
    % in results/yyyy-mm-dd/runX-THHMM-sim_case
    outdir = make_outdir(strcat("dde_sim_", run_mode)); 
else
    % placeholder empty string
    outdir = string(missing);
end

% set core parameter value set
% toy: toy values for on-off case for nondimensionalizing
% jh_opt_median_230307: Jung-Hornig simultaneous optimization performed 3/7/23
param_set = "toy";

% Status message options
% "none" = none (except what's auto-returned by matlab)
% ** "key" = some highlights (RECOMMENDED)
% "all" = all (NOT RECOMMENDED when using loops)
status_message = 'key';

% specify key times (h) for constitutive and pulse-chase
secr_end_time = 72;     % constitutive: 72h f/u for media change + treatment
pulse_length = 20/60;   % pulse = 20 minutes
chase_length = 10;      % chase = 10 hours
time_interval = 10/60;   % record results every 10 minutes
pulse_tspan = -pulse_length:time_interval:0;    % set pulse time span vector
chase_tspan = 0:time_interval:chase_length;     % set chase time span vector

%% DEFINE SYSTEM

% define species (so order of indices doesn't need to be memorized)
sp_names = ["X", "X_deg", "I", "I_deg"];

% Assign indices to species
for i = 1:length(sp_names)
    sp.(sp_names(i)) = i;
end

% DEFINE PARAMETER VALUES

switch param_set

    case "toy"

        % toy values for nondimensionalization run
        p.alpha = 1e4;
        p.beta = .1;
        p.gamma = .1;
        p.delta = .1;
        p.tau = log(2)/p.delta;

    case "jh_opt_median_230307"
        % median optimal values from Jung+Hornig optimization 3/7/23
        p.alpha = 107649;       % I production rate constant (#/cell/h)
        p.beta = .04735;       % I -> X secretion rate constant (h^-1)
        p.gamma = .1440;    % I degradation rate constant (h^-1)
        p.delta = .0204;    % X degradation rate constant (h^-1)
        p.epsilon = 0;      % X -> I internalization rate constant (h^-1)
        p.kappa = 9.005;         % transfer function variable (changes meaning when prod_decay option)
        p.tau = 1.76;     % I -> X secretion lag time (h)

end

% extract parameter names from p 
param_names = fieldnames(p);

% assign indices to parameters
for i = 1:length(param_names)
    p_ind.(param_names{i}) = i;
end

% assign number of cells - geometry taken from Hornig et al., 2000
% 5000 HUVEC cells in 0.3 mL of media
% (Jung paper didn't list number of cells? just confluency?)
% FUTURE - put these in geom_cf struct for cleaner passing/export
num_cells = 5000;   % number of cells
Vx = 3e-4;          % extracellular volume (L/well)

% Specify conversion factor: (#/cell -> pmol/well)
% (#/cell * cells/well * pmol/# = pmol/well)
% FUTURE - put these in geom_cf struct for cleaner passing/export
Nav_pmol = 6.022e11; 
conv_factor_pmol = num_cells ./ Nav_pmol;

% specify conversion factor: (pmol/well to ng/mL)
% pmol/well * mol/pmol * L/mL * g/mol * ng/g = ng/well
% ng/mL = pmol/well * mol/pmol * L/mL * g/mol * ng/g * well/L
%pmol_2_ngml = 1e-12 * 1e-3 * 1.1e5 * 1e9 / Vx;      % literature MW sFlt1: 110 kDa
% FUTURE - put these in geom_cf struct for cleaner passing/export
pmol_2_ngml = 1e-12 * 1e-3 * 9e4 * 1e9 / Vx;     % manuscript MW sFlt1: 90 kDa

% specify conversion factor: (#/cell -> ng/mL)
% FUTURE - put these in geom_cf struct for cleaner passing/export
conv_factor_ngml = conv_factor_pmol * pmol_2_ngml;

