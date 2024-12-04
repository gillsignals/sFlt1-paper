% load experimental data
load ../saved-data/jung_v2 jung_v2;     % Jung et al., 2012 uninhibited
load ../saved-data/hornig hornig;       % Hornig et al., 2000 (using basal only)
load ../saved-data/kinghorn_ctrl_summary kinghorn   % Kinghorn et al., 2023 control

% set cost options - relative cost
jung_cost_option = "omit_last_Sx_relative"; % omit 10h X point from Jung optimization
hornig_cost_option = "relative_omit3h";     % omit 3h X point from Hornig optimization
kinghorn_cost_option = "all_pts_relative";

prod_decays = ["off"];
mat_delays = ["yes"];
internalize_cases = ["no"];

n_decays = length(prod_decays);
n_mat_delays = length(mat_delays);
n_int_cases = length(internalize_cases);
num_candidate_models = n_decays * n_mat_delays * n_int_cases;

% define number of parameter sets to generate 
% (same parameters used for all candidate models)
%num_param_sets = 10; % debug
num_param_sets = 1000; % actual

num_total_runs = num_candidate_models * num_param_sets;

% initial production rate constants: log uniformly distributed
% 2/11/23: random numbers between 5e3 and 5e6
initial_alpha = 5*10.^(3*rand(num_param_sets,1) + 3);
        
% initial secretion rate constants: log uniformly distributed
% 3/7/23: random numbers between 1e-3 and 1
initial_beta = 10.^-(3*rand(num_param_sets,1));
        
% initial intracellular degradation rate constants: log uniformly distributed
% 3/7/23: random numbers between 1e-3 and 1
initial_gamma = 10.^-(3*rand(num_param_sets,1));
        
% initial intracellular degradation rate constants: log uniformly distributed
% 3/7/23: random numbers between 1e-3 and 1
initial_delta = 10.^-(3*rand(num_param_sets,1));
        
% initial secretion lag: unif dist random between 1 and 4
initial_tau = 3*rand(num_param_sets,1) + 1;
        
% create parameter matrix
param_mat = [initial_alpha, initial_beta, initial_gamma, initial_delta, initial_tau];

if save_params
    save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'pulse_tspan', ...
        'chase_tspan', 'sp', 'p', 'initial_alpha', 'initial_beta', ...
        'initial_gamma', 'initial_delta',  'initial_tau',...
         'param_mat', 'conv_factor_ngml', 'num_candidate_models', ...
        'num_param_sets', 'num_total_runs', 'jung_cost_option', ...
        'kinghorn_cost_option', 'hornig_cost_option');
end
       
% initialize counter (tracks number of bytes to update message)
%nbytes = fprintf('Testing candidate model set %d of %d...\n', 1, num_total_runs);

% initialize matrix to store optimal parameter values
optimal = zeros(num_param_sets, length(param_names));
cost = zeros(num_param_sets, 1);
exit_flag = zeros(num_param_sets, 1);

% set optimization options
lsq_options = optimoptions('lsqnonlin', 'Display', 'none');

% set optimization bounds
lb = 1e-4 * ones(length(param_names), 1);
ub = 1e7 * ones(length(param_names), 1);


