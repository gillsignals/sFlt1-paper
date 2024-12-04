% initialize array for model metadata
% columns: prod_decay, mat_delay, internalize_case, aic_param_count
model_meta = strings(num_candidate_models, 4);

% initialize output containers for simulation data
% jung_num_sim_pts: number of points in simulation output;
%               includes 2 bonus pts before pulse for plotting purposes
% vector jung_T_all(:,i) gives Jung time series for ith parameter set
% matrix Y_all(:,:,i) gives Jung simulation output matrix for ith param set
jung_num_sim_pts = length(pulse_tspan) + length(chase_tspan) + 2;
jung_T_all = ones(jung_num_sim_pts, num_total_runs) * NaN; 
jung_Y_all = ones(jung_num_sim_pts, length(fieldnames(sp)), num_total_runs) * NaN;
jung_X_frac_all = ones(jung_num_sim_pts, num_total_runs) * NaN; 
jung_I_fc_all = ones(jung_num_sim_pts, num_total_runs) * NaN; 


% jung_Si_cost_detail(i,j) = for ith parameter set, difference at jth 
%     experimental time point between simulated and observed Si metrics
jung_I_cost_detail = ones(num_total_runs, length(jung_v2.time)) * NaN;
jung_X_cost_detail = ones(num_total_runs, length(jung_v2.time)) * NaN;

% jung_Si_cost_ssd(i) = total cost from Si
%     sum of squared differences of jung_Si_cost_detail(i,:)
jung_I_cost_ssd = ones(num_total_runs, 1) * NaN;
jung_X_cost_ssd = ones(num_total_runs, 1) * NaN;

% parameters for Hornig run
hornig_secr_end_time = 72;  % time to end simulation after media change + treatment

% set time span for reporting output
hornig_opt_tspan = 0:time_interval:hornig_secr_end_time;

% initialize output containers for simulation data
% num_sim_pts: number of points in simulation output
% vector T_all(:,i) gives time series for ith parameter set
% matrix Y_all(:,:,i) gives simulation output matrix for ith param set
num_sim_pts = length(hornig_opt_tspan);
hornig_T_all = ones(num_sim_pts, num_total_runs) * NaN; 
hornig_Y_all = ones(num_sim_pts, length(fieldnames(sp)), num_total_runs) * NaN;

% hornig_Sx_cost_detail(i,j) = for ith parameter set, difference at jth 
%     experimental time point between simulated and observed Sx metrics
%     (note first 2 hornig time points omitted - treating as BLQ)
hornig_X_cost_detail = ones(num_total_runs, length(hornig.basal_time)-2) * NaN;

% hornig_Sx_cost_ssd(i) = total cost from Hornig Sx experiment
%     sum of squared differences of jung_Sx_cost_detail(i,:)
hornig_X_cost_ssd = ones(num_total_runs, 1) * NaN;

% kinghorn_X_cost_detail(i,j) = for ith parameter set, difference at jth 
%     experimental time point between simulated and observed X metrics
kinghorn_X_cost_detail = ones(num_total_runs, length(kinghorn.time)) * NaN;
kinghorn_I_cost_detail = ones(num_total_runs, length(kinghorn.time)) * NaN;

% kinghorn_X_cost_ssd(i) = total cost from Kinghorn X experiment
%     sum of squared differences of kinghorn_X_cost_detail(i,:)
kinghorn_X_cost_ssd = ones(num_total_runs, 1) * NaN;
kinghorn_I_cost_ssd = ones(num_total_runs, 1) * NaN;

% tracker of infinite sadness - which initial parameter sets don't converge?
infinite_sadness = zeros(num_total_runs,1);

% tracker of cost mismatch - lsqnonlin cost doesn't match
% results of simulations
cost_mismatch = zeros(num_total_runs,1);