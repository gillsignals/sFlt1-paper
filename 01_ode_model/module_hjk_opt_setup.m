

% load experimental data
load ../saved-data/jung_v2 jung_v2;     % Jung et al., 2012 uninhibited
load ../saved-data/hornig hornig;       % Hornig et al., 2000 (using basal only)
load ../saved-data/kinghorn_ctrl_summary kinghorn   % Kinghorn et al., 2023 control

% set cost options - relative cost
jung_cost_option = "omit_last_Sx_relative"; % omit 10h X point from Jung optimization
hornig_cost_option = "relative_omit3h";     % omit 3h X point from Hornig optimization
kinghorn_cost_option = "all_pts_relative";

% define number of parameter sets to generate 
% (same parameters used for all candidate models)
switch run_mode
    case "base_constitutive"
        num_param_sets = 1;
    case "base_pulse_chase"
        num_param_sets = 1;
    case "opt_1"
        num_param_sets = 1;
    case "opt_100"
        %num_param_sets = 10; % debug
        num_param_sets = 100; % actual
end

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
        
% create parameter matrix
param_mat = [initial_alpha, initial_beta, initial_gamma, initial_delta];

if save_params
    save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'pulse_tspan', ...
        'chase_tspan', 'sp', 'p', 'initial_alpha', 'initial_beta', 'initial_gamma', ...
        'initial_delta', 'param_mat', 'conv_factor_ngml', 'num_param_sets', ...
        'jung_cost_option', 'kinghorn_cost_option', 'hornig_cost_option');
end

% set optimization options
lsq_options = optimoptions('lsqnonlin', 'Display', 'none');

% set optimization bounds
lb = [1e3, 1e-6, 1e-6, 1e-6];
ub = [1e7, 1, 1, 1];

% initialize matrices to store optimal parameter values
optimal = zeros(num_param_sets, length(param_names));
cost = zeros(num_param_sets, 1);
exit_flag = zeros(num_param_sets, 1);

% initialize output containers for simulation data
% jung_num_sim_pts: number of points in simulation output;
%               includes 2 bonus pts before pulse for plotting purposes
% vector jung_T_all(:,i) gives Jung time series for ith parameter set
% matrix Y_all(:,:,i) gives Jung simulation output matrix for ith param set
jung_num_sim_pts = length(pulse_tspan) + length(chase_tspan) + 2;
jung_T_all = ones(jung_num_sim_pts, num_param_sets) * NaN; 
jung_Y_all = ones(jung_num_sim_pts, length(fieldnames(sp)), num_param_sets) * NaN;
jung_X_frac_all = ones(jung_num_sim_pts, num_param_sets) * NaN; 
jung_I_fc_all = ones(jung_num_sim_pts, num_param_sets) * NaN; 


% jung_Si_cost_detail(i,j) = for ith parameter set, difference at jth 
%     experimental time point between simulated and observed Si metrics
jung_I_cost_detail = ones(num_param_sets, length(jung_v2.time)) * NaN;
jung_X_cost_detail = ones(num_param_sets, length(jung_v2.time)) * NaN;

% jung_Si_cost_ssd(i) = total cost from Si
%     sum of squared differences of jung_Si_cost_detail(i,:)
jung_I_cost_ssd = ones(num_param_sets, 1) * NaN;
jung_X_cost_ssd = ones(num_param_sets, 1) * NaN;

% set time span for reporting output
hornig_opt_tspan = 0:time_interval:secr_end_time;

% initialize output containers for simulation data
% num_sim_pts: number of points in simulation output
% vector T_all(:,i) gives time series for ith parameter set
% matrix Y_all(:,:,i) gives simulation output matrix for ith param set
num_sim_pts = length(hornig_opt_tspan);
hornig_T_all = ones(num_sim_pts, num_param_sets) * NaN; 
hornig_Y_all = ones(num_sim_pts, length(fieldnames(sp)), num_param_sets) * NaN;

% hornig_Sx_cost_detail(i,j) = for ith parameter set, difference at jth 
%     experimental time point between simulated and observed Sx metrics
%     (note first 2 hornig time points omitted - treating as BLQ)
hornig_X_cost_detail = ones(num_param_sets, length(hornig.basal_time)-2) * NaN;

% hornig_Sx_cost_ssd(i) = total cost from Hornig Sx experiment
%     sum of squared differences of jung_Sx_cost_detail(i,:)
hornig_X_cost_ssd = ones(num_param_sets, 1) * NaN;

% kinghorn_X_cost_detail(i,j) = for ith parameter set, difference at jth 
%     experimental time point between simulated and observed X metrics
kinghorn_X_cost_detail = ones(num_param_sets, length(kinghorn.time)) * NaN;
kinghorn_I_cost_detail = ones(num_param_sets, length(kinghorn.time)) * NaN;

% kinghorn_X_cost_ssd(i) = total cost from Kinghorn X experiment
%     sum of squared differences of kinghorn_X_cost_detail(i,:)
kinghorn_X_cost_ssd = ones(num_param_sets, 1) * NaN;
kinghorn_I_cost_ssd = ones(num_param_sets, 1) * NaN;

% tracker of infinite sadness - which initial parameter sets don't converge?
infinite_sadness = zeros(num_param_sets,1);

% tracker of cost mismatch - lsqnonlin cost doesn't match
% results of simulations
cost_mismatch = zeros(num_param_sets,1);