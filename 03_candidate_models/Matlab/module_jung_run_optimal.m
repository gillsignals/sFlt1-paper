% simulate Jung results with optimal parameter set to store time course output
            [jung_sim_time, jung_sim_y, jung_sim_X_frac, jung_sim_I_fc] = ...
                sim_pulse_chase_dde(prod_decay, mat_delay, internalize_case,  ...
                sp, p, pulse_length, chase_length, time_interval, conv_factor_ngml);
                 
            % store time courses in a single data structure
            jung_T_all(:,out_index) = jung_sim_time;
            jung_Y_all(:,:,out_index) = jung_sim_y';
            jung_X_frac_all(:,out_index) = jung_sim_X_frac;
            jung_I_fc_all(:,out_index) = jung_sim_I_fc;

            % calculate cost of Jung fit
            for j = 1:length(jung_v2.time)
        
                % find simulated time point closest to experimental time point
                % [~, tindex] = min(abs(T_chase - T_exp(j))); EQUIVALENT TO NEXT LINE
                tindex = dsearchn(jung_sim_time', jung_v2.time(j));
        
                % calculate costs depending on stated cost_option case
                switch jung_cost_option
        
                    % (simulated - experimental) for all Sx, Si fold changes
                    case 'all_pts_absolute'
                        jung_X_cost_detail(out_index,j) = jung_sim_X_frac(tindex) - jung_v2.media(j);
                        jung_I_cost_detail(out_index,j) = jung_sim_I_fc(tindex) - jung_v2.lysate(j);
                               
                    % (simulated - experimental)/experimental for all Sx, Si fold changes
                    case 'all_pts_relative'
                        if jung_v2.media(j) ~=0
                            jung_X_cost_detail(out_index,j) = (jung_sim_X_frac(tindex) - ...
                                jung_v2.media(j))/jung_v2.media(j);
                        end
                        jung_I_cost_detail(out_index,j) = (jung_sim_I_fc(tindex) - ...
                            jung_v2.lysate(j))/jung_v2.lysate(j);
        
                    % (simulated - experimental) omitting last Sx data point
                    case 'omit_last_Sx_absolute'
                        if jung_v2.time(j) ~= 10
                            % error normalized to observed value
                            jung_X_cost_detail(out_index,j) = jung_sim_X_frac(tindex) - jung_v2.media(j);
                        end
                        jung_I_cost_detail(out_index,j) = jung_sim_I_fc(tindex) - jung_v2.lysate(j);
        
                    % (simulated - experimental)/experimental for all Sx, Si fold changes
                    case 'omit_last_Sx_relative'
                        if jung_v2.time(j) ~= 10 & jung_v2.media(j) ~= 0
                            % error normalized to observed value
                            jung_X_cost_detail(out_index,j) = (jung_sim_X_frac(tindex) - jung_v2.media(j))/jung_v2.media(j);
                        end
                        jung_I_cost_detail(out_index,j) = (jung_sim_I_fc(tindex) - jung_v2.lysate(j))/jung_v2.lysate(j);
                end
        
            end      