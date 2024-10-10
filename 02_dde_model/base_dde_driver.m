%% Analysis of candidate DDE models of sFLT1 secretion optimized to Hornig, Jung and Kinghorn data

clear;
close all;

% set run mode
% "base_on_off": run one simulation with equal periods of production and none
% "vary_tau_on_off": run simulations with several values of tau
% "base_model_opt": run optimization to experimental data 1000x
run_mode = "vary_tau_on_off";

% Announce which case will run
disp("Running case " + run_mode + "...")

% Data saving options
save_params = 0;    % Option to save input parameters; 1 = yes, 0 = no
save_data = 0;  % Option to save data; 1 = yes, 0 = no
save_figs = 0;  % Option to save figures; 1 = yes, 0 = no

module_setup;

%% CASES
switch run_mode

    case 'base_on_off'

        % characteristic time unit: I half-time to steady state
        t50_X = log(2)/p.delta;

        % characteristic I unit: I at steady state
        I_ss = p.alpha / (p.beta + p.gamma);

        % characteristic X unit: X at steady state
        X_ss = p.beta * I_ss / p.delta;

        % only 1 parameter set
        num_param_sets = 1;

        % lengthen pulse and chase to see long-term dynamics
        num_time_units = 12; 
        pulse_length = num_time_units * t50_X;      % set pulse based on X half-life
        chase_length = num_time_units * t50_X;      % set chase based on X half-life
        time_interval = t50_X / 10;         % record 10 results per time unit
        pulse_tspan = -pulse_length:time_interval:0;    % set pulse time span vector
        chase_tspan = 0:time_interval:chase_length;     % set chase time span vector

        % save simulation input parameters before running
        if save_params
            save(fullfile(outdir, 'inputs', 'on_off_input.mat'), 'num_cells', 'pulse_length', ...
                'chase_length', 'time_interval', 'sp', 'p', 'conv_factor_ngml', 'num_time_units');
        end

        % run simulation
        [pc_sim_time, pc_sim_Y, pc_sim_X_norm, pc_sim_I_norm] = sim_on_off_dde(sp, ...
            p, pulse_length, chase_length, time_interval, t50_X, I_ss, X_ss, conv_factor_ngml);
    
        % save time course, parameter values, options
        if save_data
            save(fullfile(outdir, 'results', 'on_off_sim.mat'), ...
                'pc_sim_time', 'pc_sim_Y', 'pc_sim_X_norm', 'pc_sim_I_norm', ...
                'pulse_length', 'chase_length', 'time_interval');
        end

        % initialize plot - absolute amounts
        abs_fig = figure("Name","On-off time course plot - absolute amounts");

        % Extracellular sFlt1 (X) absolute amount
        subplot(2,1,1);
        plot(pc_sim_time, pc_sim_X_norm, 'LineWidth', 3)
        xlim([-14, 12])
        xline(-12, '-k', "Switch on", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Switch off", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Extracellular sFLT1 Amount (X)")
        xlabel("Time (h)")
        ylabel("X (normalized)")
        
        % Intracellular sFlt1 (I) absolute amount
        subplot(2,1,2)
        plot(pc_sim_time, pc_sim_I_norm, 'LineWidth', 3)
        xlim([-14, 12])
        xline(-12, '-k', "Switch on", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Switch off", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Intracellular sFLT1 Amount (I)")
        xlabel("Time (h)")
        ylabel("I (normalized)")


    case 'vary_tau_on_off'

        tic

        % characteristic time unit: I half-time to steady state
        t50_X = log(2)/p.delta;

        % specify tau values to test
        taus = ([1e-6,.25:0.25:2]) * t50_X;
        num_param_sets = length(taus);

        % set pulse and chase lengths to see long-term dynamics
        num_time_units = 12; 
        pulse_length = num_time_units * t50_X;      % set pulse based on X half-life
        chase_length = num_time_units * t50_X;      % set chase based on X half-life
        time_interval = t50_X / 10;         % record 10 results per time unit
        pulse_tspan = -pulse_length:time_interval:0;    % set pulse time span vector
        chase_tspan = 0:time_interval:chase_length;     % set chase time span vector
        num_sim_pts = length(pulse_tspan) + length(chase_tspan) + 2;
        
        % save simulation input parameters before running
        if save_params
            save(fullfile(outdir, 'inputs', 'on_off_taus_input.mat'), 'num_cells',  ...
                'pulse_length', 'chase_length', 'time_interval', 'sp', 'p', ...
                'conv_factor_ngml', 'num_time_units', 'taus');
        end

        % initialize output arrays and vectors
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), num_param_sets) * NaN;
        X_norm_all = ones(num_sim_pts, num_param_sets) * NaN;
        I_norm_all = ones(num_sim_pts, num_param_sets) * NaN;

        % characteristic I unit: I at steady state
        I_ss = p.alpha / (p.beta + p.gamma);

        % characteristic X unit: X at steady state
        X_ss = p.beta * I_ss / p.delta;

        % run for each tau value
        for i = 1:length(taus)

            p.tau = taus(i);

            [pc_sim_time, pc_sim_Y, pc_sim_X_norm, pc_sim_I_norm] = sim_on_off_dde(sp, ...
                p, pulse_length, chase_length, time_interval, t50_X, I_ss, X_ss, conv_factor_ngml);
        
            Y_all(:,:,i) = pc_sim_Y';
            X_norm_all(:,i) = pc_sim_X_norm;
            I_norm_all(:,i) = pc_sim_I_norm;

        end

        % save time course, parameter values, options
        if save_data
            save(fullfile(outdir, 'results', 'on_off_sim_taus.mat'), ...
                'pc_sim_time', 'Y_all', 'X_norm_all', 'I_norm_all', ...
                'pulse_length', 'chase_length', 'time_interval');
        end

        toc

    case 'base_model_opt'

        tic

        if ismember(status_message, [{'all'}])
            disp("Running case " + run_mode + "...")
        end   

        module_hjk_opt_setup;

        module_init_hjk_arrays;

        % generate candidate models and run n=num_param_sets simulations for each 
        for a = 1:n_decays
            for b = 1:n_mat_delays
                for c = 1:n_int_cases

                    % set flux term options for candidate model equations
                    prod_decay = prod_decays(a);
                    mat_delay = mat_delays(b);
                    internalize_case = internalize_cases(c);

                    % calculate parameters / degrees of freedom: 
                    % 4 base params + up to 3 additional + 1 for variance
                    aic_param_count = 5 + (prod_decay ~= "off") + ...
                        (mat_delay == "yes") + (internalize_case == "yes");

                    model_num = (a-1) * n_mat_delays * n_int_cases + ...
                        (b-1) * n_int_cases + c;

                    model_meta(model_num, :) = [prod_decay, mat_delay, ...
                        internalize_case, aic_param_count];

                    % optimize candidate model starting from each parameter set
                    for i = 1:num_param_sets

                        % calculate index of this model run
                        out_index = (a-1) * n_mat_delays * n_int_cases * num_param_sets + ...
                            (b-1) * n_int_cases * num_param_sets + ...
                            (c-1) * num_param_sets + i;
    
                        % update counter
                        %fprintf(repmat('\b',1,nbytes))
                        nbytes = fprintf('Optimizing model %d of %d...\n', ...
                            out_index, num_candidate_models * num_param_sets);
    
                        % set parameters to the ith set
                        k_init = param_mat(i,:);
                        for j = 1:length(param_names)
                           p.(param_names{j}) = k_init(j);
                        end

                        %% optimize

            % define anonymous cost function to pass additional arguments
            cost_fxn = @(k_init) hjk_combo_cost(k_init, kinghorn.time, kinghorn.X_frac_24h, ...            
                kinghorn.I_fc_0h, jung_v2.time, jung_v2.media, jung_v2.lysate, ...
                hornig.basal_time, hornig.basal_Sx_ng_ml, prod_decay, mat_delay, ...
                internalize_case, sp, p, param_names, conv_factor_ngml,  pulse_length, ...
                chase_length, hornig_secr_end_time, time_interval, jung_cost_option, ...
                hornig_cost_option, kinghorn_cost_option, outdir, status_message);

            % run nonlinear least squares optimization
            try
                [optimal(out_index,:), cost(out_index), ~, exit_flag(out_index)] = ...
                    lsqnonlin(cost_fxn, k_init, lb, ub, lsq_options);
            catch
                infinite_sadness(out_index) = 1;
                fprintf("Failed parameter sets: %d \n", sum(infinite_sadness));
                continue
            end

            % set parameters to optimized values
            for j = 1:length(param_names)
                p.(param_names{j}) = optimal(out_index,j);
            end

            try
                module_jung_run_optimal;
                module_constitutive_run_optimal;
            catch
                continue;
            end

            % calculate total cost of fit to each dataset
            jung_I_cost_ssd(out_index) = sum(jung_I_cost_detail(out_index,:).^2, "omitnan");
            jung_X_cost_ssd(out_index) = sum(jung_X_cost_detail(out_index,:).^2, "omitnan");    
            hornig_X_cost_ssd(out_index) = sum(hornig_X_cost_detail(out_index,:).^2, "omitnan");
            kinghorn_X_cost_ssd(out_index) = ...
                sum(kinghorn_X_cost_detail(out_index,:).^2, "omitnan");
            kinghorn_I_cost_ssd(out_index) = ...
                sum(kinghorn_I_cost_detail(out_index,:).^2, "omitnan");

            % test that cost from lsqnonlin matches calculated cost of fit
            % (error of less than cost/1000)
            cost_error = cost(out_index) - jung_I_cost_ssd(out_index) - ...
                    jung_X_cost_ssd(out_index) - hornig_X_cost_ssd(out_index) - ...
                    kinghorn_X_cost_ssd(out_index) - kinghorn_I_cost_ssd(out_index);
            if cost_error / cost(out_index) > .001
                warning("Cost from lsqnonlin does not match cost observed with " + ...
                    "optimization run %d. \n", out_index);

                cost_mismatch(out_index) = 1;
            end

                    end
                end
            end
        end

disp("Saving optimized parameters and time courses, if desired...")

        % save optimal parameters and time course data for this cost/decay option set
        if save_data
            save(fullfile(outdir, 'results', 'hjk_dde_opt_1000.mat'), ...
                'optimal', 'cost', 'exit_flag', 'infinite_sadness', 'cost_mismatch', ...
                'jung_T_all', 'jung_Y_all', 'jung_X_frac_all', 'jung_I_fc_all', ...
                'jung_X_cost_detail', 'jung_X_cost_ssd', 'jung_I_cost_detail', ...
                'jung_I_cost_ssd', 'hornig_T_all', 'hornig_Y_all', ...
                'hornig_X_cost_detail', 'hornig_X_cost_ssd', 'model_meta', ...
                'kinghorn_X_cost_ssd', 'kinghorn_X_cost_detail', ...
                'kinghorn_I_cost_ssd', 'kinghorn_I_cost_detail');

            writematrix(model_meta, fullfile(outdir, 'inputs', 'hjk_dde_opt_meta.csv'))
        end

        % note hornig_Y_all becomes too large to save in one file 
        % at large num_parameter_sets - save in chunks if needed
        %   hornig_spX_all = hornig_Y_all(:, sp.X, :);
        %   hornig_spI_all = hornig_Y_all(:, sp.I, :);
        %   save(fullfile(outdir, 'results', 'jh_opt_candidates_1000_hornig_sp.mat'), 'hornig_spX_all', 'hornig_spI_all')

        
        toc

    
        
end