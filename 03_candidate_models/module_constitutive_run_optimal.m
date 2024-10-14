%% RUN CONSTITUTIVE SIMULATION WITH OPTIMAL PARAMETER SET

% simulate constitutive results with optimal parameter set to store time course output
[constitutive_sim_time, constitutive_sim_Y] = sim_secr_dde_v2(sp, p, "on", mat_delay, ...
    internalize_case, hornig_secr_end_time, time_interval, conv_factor_ngml, status_message);

% store time courses in a single data structure
hornig_T_all(:,out_index) = constitutive_sim_time;
hornig_Y_all(:,:,out_index) = constitutive_sim_Y'; 

%% CALCULATE HORNIG COST

% calculate cost of Hornig fit
% generate vector of Hornig costs
for j=3:length(hornig.basal_time)   % removing early concs potentially below detection limit
            
    % find simulated time point closest to experimental time point
    % [~, tindex] = min(abs(T_chase - T_exp(j))); EQUIVALENT TO NEXT LINE
    t_index = dsearchn(constitutive_sim_time', hornig.basal_time(j));
            
    % calculate costs depending on stated cost_option case
    switch hornig_cost_option
                
        case 'absolute_omit3h'
            % calculate difference between simulated and observed Sx
            hornig_X_cost_detail(out_index,j-2) = constitutive_sim_Y(sp.X, t_index) - ...
                hornig.basal_Sx_ng_ml(j);
                
        case 'relative_omit3h'
            % calculate normalized difference between simulated and observed Sx
            hornig_X_cost_detail(out_index,j-2) = (constitutive_sim_Y(sp.X, t_index) - ...
                hornig.basal_Sx_ng_ml(j)) ./ hornig.basal_Sx_ng_ml(j);
    end
            
end

%% CALCULATE KINGHORN COST

% convert X to fraction of 24h signal
sim_X_24h = constitutive_sim_Y(sp.X, dsearchn(constitutive_sim_time', 24));
sim_X_frac_24h = constitutive_sim_Y(sp.X,:) ./ sim_X_24h;

% convert I to fold change versus 0h
sim_I_0h = constitutive_sim_Y(sp.I, dsearchn(constitutive_sim_time', 0));
sim_I_fc_0h = constitutive_sim_Y(sp.I,:) ./ sim_I_0h;
    
% generate vectors of Kinghorn costs for X and I
for j = 1:length(kinghorn.time)

    % find simulated time point closest to experimental time point
    % [~, tindex] = min(abs(T_chase - T_exp(j))); EQUIVALENT TO NEXT LINE
    tindex = dsearchn(constitutive_sim_time', kinghorn.time(j));

    % calculate costs depending on stated cost_option case
    switch kinghorn_cost_option

        % (simulated - experimental) for all X, I fold changes
        case 'all_pts_absolute'
            kinghorn_X_cost_detail(j) = sim_X_frac_24h(tindex) - kinghorn.X_frac_24h(j);
            kinghorn_I_cost_detail(j) = sim_I_fc_0h(tindex) - kinghorn.I_fc_0h(j);
                       
        % (simulated - experimental)/experimental for all X, I fold changes
        case 'all_pts_relative'
            if kinghorn.X_frac_24h(j) ~= 0 & sim_X_frac_24h(tindex) ~= 0
                kinghorn_X_cost_detail(out_index, j) = (sim_X_frac_24h(tindex) - ...
                    kinghorn.X_frac_24h(j))/kinghorn.X_frac_24h(j);
            end
            kinghorn_I_cost_detail(out_index, j) = (sim_I_fc_0h(tindex) - ...
                kinghorn.I_fc_0h(j))/kinghorn.I_fc_0h(j);
    end
end
    
