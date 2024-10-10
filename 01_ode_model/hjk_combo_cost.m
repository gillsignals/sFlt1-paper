function cost = hjk_combo_cost(k_try, kinghorn_T_exp, kinghorn_X_exp, kinghorn_I_exp, ...
    jung_T_exp, jung_X_exp, jung_I_exp, hornig_T_exp, hornig_X_exp, sp, p, param_names, ...
    conv_factor, pulse_length, chase_length, secr_end_time, time_interval, ...
    jung_cost_option, hornig_cost_option, kinghorn_cost_option)


%% KEY QUESTION - how to balance sets of costs?
% Jung has both X, I, 5 time points each (omitting trivial 0h)
% (may need to count as 4 each because normalized to 8h, so that's also trivial)
% Hornig has 6 nonzero X (omitting 0h and 3h)
% Kinghorn has 40 I points (omitting trivial 0h) and 40 X points (omitting trivial 24h)
% Or, reduce to 6 each by averaging across each time
% how to weight these?
%   different scales -> relative sensitivity
%   different numbers of data points -> after talking to Feilim, first pass
%       will not use any correction


% hjk_combo_cost  Cost function to optimize trafficking parameters to 
%               simultaneously fit Jung et al., 2012 pulse-chase data, 
%               one Hornig et al., 2000 dataset, and control data from
%               Kinghorn et al., 2023, within lsqnonlin using specified equations
%               
%
% EXAMPLE USAGE (within lsqnonlin)
% --------------------------------
%
% % define anonymous cost function to pass additional arguments
% cost_fxn = @(k_init) jung_hornig_combo_cost(k_init, kinghorn.time,
%       kinghorn.X_frac_24h, kinghorn.I_fc_0h, jung_v2.time, jung_v2.media, ...
%       jung_v2.lysate, hornig.basal_time, hornig.basal_Sx_ng_ml, sp, p,
%       param_names, conv_factor_ngml, pulse_length, chase_length, 
%       hornig_secr_end_time, time_interval, cost_option, status_message);
%
%
% % run nonlinear least squares optimization
% % (lb = lower bounds, ub = upper bounds, lsq_options = lsqnonlin options)
% optimal(x, :) = lsqnonlin(cost_fxn, k_init, lb, ub, lsq_options);
%
% OUTPUT (to lsqnonlin)
% ---------------------
% cost:         Current cost: sum of squared difference between
%               experimental and simulated values across all data points
%
% INPUT
% -----
% k_try:        Vector of initial/current guesses for parameter values
% jung_T_exp:   Vector of experimental times from Jung pulse-chase
% jung_X_exp:  Vector of experimental extracellular fold changes in Jung
% jung_I_exp:  Vector of experimental intracellular fold changes in Jung
% hornig_T_exp: Vector of experimental times from Hornig expt
% hornig_X_exp:    Vector of Hornig expt extracellular sFlt1 (ng/mL)
% prod_case:   String noting type of production decay to use in chase
% sp:           Struct mapping species names to index values
% p:            Struct holding parameter values; updated throughout opt
% param_names:  Vector containing names of struct p, for mapping k_try to p
% conv_factor:  Conversion factor for #/cell ->ng/ml
% pulse_length: Scalar representing duration of pulse (h)
% chase_tspan:  Scalar representing duration of chase (h)
% hornig_secr_end_time: Duration of Hornig constitutive secretion
% time_interval: Scalar representing how often to return simulated data points (h)
% ode_options:  ODE solver options (previously specified via odeset)
% cost_option:  specifies which points to include/omit and whether to scale
% eqn_case:     String specifying which set of equations to use (see dde_driver)
% outdir: directory prefix to be used for saving results
% status_message: Controls messages returned by Matlab (see dde_driver)

    %% SET PARAMETERS

    % package current parameter estimates into p struct
    for i = 1:length(param_names)
        p.(param_names{i}) = k_try(i);
    end


try
    
    %% JUNG SIMULATION AND COST
    
    % initialize cost vectors to 0, removing last point as potential outlier
    jung_X_cost = zeros(length(jung_T_exp) - 1,1);
    jung_I_cost = zeros(length(jung_T_exp) - 1,1);

    % run simulation - DO NOT save parameters
    [jung_sim_time, ~, sim_X_frac, sim_I_fc] = sim_pulse_chase_ode(sp, p, ...
    pulse_length, chase_length, time_interval, conv_factor);

    % generate vectors of Jung costs for X and I
    for j = 1:length(jung_T_exp)

        % find simulated time point closest to experimental time point
        % [~, tindex] = min(abs(T_chase - T_exp(j))); EQUIVALENT TO NEXT LINE
        tindex = dsearchn(jung_sim_time', jung_T_exp(j));

        % calculate costs depending on stated cost_option case
        switch jung_cost_option

            % (simulated - experimental) for all X, I fold changes
            case 'all_pts_absolute'
                jung_X_cost(j) = sim_X_frac(tindex) - jung_X_exp(j);
                jung_I_cost(j) = sim_I_fc(tindex) - jung_I_exp(j);
                       
            % (simulated - experimental)/experimental for all X, I fold changes
            case 'all_pts_relative'
                if jung_X_exp(j) ~=0
                    jung_X_cost(j) = (sim_X_frac(tindex) - jung_X_exp(j))/jung_X_exp(j);
                end
                jung_I_cost(j) = (sim_I_fc(tindex) - jung_I_exp(j))/jung_I_exp(j);

            % (simulated - experimental) omitting last Sx data point
            case 'omit_last_Sx_absolute'
                if jung_T_exp(j) ~= 10
                    % error normalized to observed value
                    jung_X_cost(j) = sim_X_frac(tindex) - jung_X_exp(j);
                end
                jung_I_cost(j) = sim_I_fc(tindex) - jung_I_exp(j);

            % (simulated - experimental)/experimental for all X, I fold changes
            case 'omit_last_Sx_relative'
                if jung_T_exp(j) ~= 10 & jung_X_exp(j) ~= 0
                    % error normalized to observed value
                    jung_X_cost(j) = (sim_X_frac(tindex) - jung_X_exp(j))/jung_X_exp(j);
                end
                jung_I_cost(j) = (sim_I_fc(tindex) - jung_I_exp(j))/jung_I_exp(j);
        end

    end

    %% CONSTITUTIVE SIMULATION AND HORNIG COST
    
    % initialize cost vector to 0
    hornig_X_cost = zeros(length(hornig_T_exp) - 2,1); % removing first 2 points bc potentially BLQ
    kinghorn_X_cost = zeros(length(kinghorn_T_exp),1);
    kinghorn_I_cost = zeros(length(kinghorn_T_exp),1);

        % run simulation to steady state followed by media change and 72h followup
        [constitutive_sim_time, constitutive_sim_Y] = sim_secr_ode(sp, p, ...
            secr_end_time, time_interval, conv_factor);
        
        % generate vector of Hornig costs
        for j=3:length(hornig_T_exp)   % removing early concs potentially below detection limit
        
            % find simulated time point closest to experimental time point
            % [~, tindex] = min(abs(T_chase - T_exp(j))); EQUIVALENT TO NEXT LINE
            t_index = dsearchn(constitutive_sim_time', hornig_T_exp(j));
        
            % calculate costs depending on stated cost_option case
            switch hornig_cost_option
            
                case 'absolute_omit3h'
                    % calculate difference between simulated and observed Sx
                    hornig_X_cost(j-2) = constitutive_sim_Y(t_index, sp.X) - hornig_X_exp(j);
            
                case 'relative_omit3h'
                    % calculate normalized difference between simulated and observed Sx
                    hornig_X_cost(j-2) = (constitutive_sim_Y(t_index, sp.X) - hornig_X_exp(j)) ...
                        ./ hornig_X_exp(j);
            end
        
        end

        %% KINGHORN COST (with Hornig simulation)

        % convert X to fraction of 24h signal
    sim_X_24h = constitutive_sim_Y(dsearchn(constitutive_sim_time', 24), sp.X);
    sim_X_frac_24h = constitutive_sim_Y(:,sp.X) ./ sim_X_24h;

    % convert I to fold change versus 0h
    sim_I_0h = constitutive_sim_Y(dsearchn(constitutive_sim_time', 0), sp.I);
    sim_I_fc_0h = constitutive_sim_Y(:,sp.I) ./ sim_I_0h;
    
    
    % generate vectors of Kinghorn costs for X and I
    for j = 1:length(kinghorn_T_exp)

        % find simulated time point closest to experimental time point
        % [~, tindex] = min(abs(T_chase - T_exp(j))); EQUIVALENT TO NEXT LINE
        tindex = dsearchn(constitutive_sim_time', kinghorn_T_exp(j));

        % calculate costs depending on stated cost_option case
        switch kinghorn_cost_option

            % (simulated - experimental) for all X, I fold changes
            case 'all_pts_absolute'
                kinghorn_X_cost(j) = sim_X_frac_24h(tindex) - kinghorn_X_exp(j);
                kinghorn_I_cost(j) = sim_I_fc_0h(tindex) - kinghorn_I_exp(j);
                       
            % (simulated - experimental)/experimental for all X, I fold changes
            case 'all_pts_relative'
                if sim_X_frac_24h(tindex) ~=0
                    kinghorn_X_cost(j) = (sim_X_frac_24h(tindex) - ...
                        kinghorn_X_exp(j))/kinghorn_X_exp(j);
                end
                kinghorn_I_cost(j) = (sim_I_fc_0h(tindex) - ...
                    kinghorn_I_exp(j))/kinghorn_I_exp(j);
        end
    end


    %% COMBINE COSTS

    % combine X costs and I costs - minimize in same lsq call
    cost = [jung_X_cost', jung_I_cost', hornig_X_cost', kinghorn_X_cost', kinghorn_I_cost'];
    

catch
    cost = NaN;
    

end
