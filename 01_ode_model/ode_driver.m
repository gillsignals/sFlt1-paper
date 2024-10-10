%% Analysis of ODE models of sFLT1 secretion optimized to Hornig, Jung and Kinghorn data

clear;
close all;

% set run mode
% "base_constitutive": run one constitutive secretion simulation with
%       default parameter values
% "base_on_off": run one simulation with equal periods of production and none
% "base_pulse_chase": run one pulse-chase simulation with default parameter values
% "extended_pulse_chase": run one pulse-chase simulation for longer time
% "opt_1": run one optimization starting from default parameter values
% "opt_100": run 100 optimizations from randomly sampled initial values
run_mode = "opt_100";

% Announce which case will run
disp("Running case " + run_mode + "...")

% Data saving options
save_params = 0;    % Option to save input parameters; 1 = yes, 0 = no
save_data = 0;  % Option to save data; 1 = yes, 0 = no
save_figs = 0;  % Option to save figures; 1 = yes, 0 = no

% make output directory, assign time intervals, and
%   specify default parameter values (kinetic, geometric, conversion factors)
module_setup;

%% CASES
switch run_mode

    case 'base_constitutive'

        module_hjk_opt_setup;

        % save simulation input parameters before running
        if save_params
            save(fullfile(outdir, 'inputs', 'constitutive_input.mat'), 'num_cells', ...
                'time_interval', 'sp', 'p', 'conv_factor_ngml');
        end
        
        % run simulation to steady state followed by media change and 72h followup
        [const_sim_time, const_sim_Y, const_sim_X_frac, const_sim_I_fc] = ...
            sim_secr_ode(sp, p, secr_end_time, time_interval, conv_factor_ngml);

        % save time course, parameter values, options
        if save_data
            save(fullfile(outdir, 'results', 'const_sim.mat'), ...
                'const_sim_time', 'const_sim_Y', 'const_sim_X_frac', 'const_sim_I_fc', ...
                'secr_end_time', 'time_interval', 'conv_factor_ngml');
        end

        % Extracellular sFlt1 (Sx) amount, ng/mL
        fig = figure;
        subplot(2,1,1);
        plot(const_sim_time, const_sim_Y(:,sp.X), 'LineWidth', 3)
        hold on;
        scatter(hornig.basal_time, hornig.basal_Sx_ng_ml)
        hold off;
        title("Extracellular sFLT1 Amount (X)")
        xlabel("Time (h)")
        ylabel("Amount (ng/ml)")
        
        % Intracellular sFlt1 (Si) amount
        subplot(2,1,2)
        plot(const_sim_time, const_sim_Y(:,sp.I), 'LineWidth', 3)
        title("Intracellular sFLT1 Amount (I)")
        xlabel("Time (h)")
        ylabel("Amount (#/cell)")

        % save figure - absolute amounts of species
        if save_figs
            saveas(fig, fullfile(outdir, 'figures', ...
                strcat(run_mode, '.png')))
        end

    case 'base_on_off'

        % set parameter values
        p.alpha = 1e4;
        p.beta = .1;
        p.gamma = .1;
        p.delta = .1;

        % characteristic time unit: X half-time to steady state
        t50_X = log(2)/p.delta;

        % characteristic I unit: I at steady state
        I_ss = p.alpha / (p.beta + p.gamma);

        % characteristic X unit: X at steady state
        X_ss = p.beta * I_ss / p.delta;

        % only 1 parameter set
        num_param_sets = 1;

        % lengthen pulse and chase to see long-term dynamics
        num_time_units = 6; 
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

        [pc_sim_time, pc_sim_Y, pc_sim_X_norm, pc_sim_I_norm] = sim_on_off_ode(sp, ...
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
        xlim([-2-num_time_units, num_time_units])
        xline(-12, '-k', "Switch on", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Switch off", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Extracellular sFLT1 Amount (X)")
        xlabel("Time (h)")
        ylabel("X (normalized)")
        
        % Intracellular sFlt1 (I) absolute amount
        subplot(2,1,2)
        plot(pc_sim_time, pc_sim_I_norm, 'LineWidth', 3)
        xlim([-2-num_time_units, num_time_units])
        xline(-12, '-k', "Switch on", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Switch off", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Intracellular sFLT1 Amount (I)")
        xlabel("Time (h)")
        ylabel("I (normalized)")

    case 'base_pulse_chase'

        module_hjk_opt_setup;

        % save simulation input parameters before running
        if save_params
            save(fullfile(outdir, 'inputs', 'pulse_chase_input.mat'), 'num_cells', 'pulse_length', ...
                'chase_length', 'time_interval', 'sp', 'p', 'conv_factor_ngml');
        end

        [pc_sim_time, pc_sim_Y, pc_sim_X_frac, pc_sim_I_fc] = sim_pulse_chase_ode(sp, ...
            p, pulse_length, chase_length, time_interval, conv_factor_ngml);
    
        % save time course, parameter values, options
        if save_data
            save(fullfile(outdir, 'results', 'pc_sim.mat'), ...
                'pc_sim_time', 'pc_sim_Y', 'pc_sim_X_frac', 'pc_sim_I_fc', ...
                'pulse_length', 'chase_length', 'time_interval');
        end

        % initialize plot - absolute amounts
        abs_fig = figure("Name","Pulse-chase time course plot - absolute amounts");

        % Extracellular sFlt1 (X) absolute amount
        subplot(2,1,1);
        plot(pc_sim_time, pc_sim_Y(:,sp.X), 'LineWidth', 3)
        xlim([-3*pulse_length, chase_length])
        xline(-pulse_length, '-k', "Start pulse", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Start chase", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Extracellular sFLT1 Amount (X)")
        xlabel("Time (h)")
        ylabel("Amount of X (ng/mL)")
        
        % Intracellular sFlt1 (I) absolute amount
        subplot(2,1,2)
        plot(pc_sim_time, pc_sim_Y(:,sp.I), 'LineWidth', 3)
        xlim([-3*pulse_length, chase_length])
        xline(-pulse_length, '-k', "Start pulse", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Start chase", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Intracellular sFLT1 Amount (I)")
        xlabel("Time (h)")
        ylabel("Amount of I (#/cell)")

        % save figure - absolute amounts of species
        if save_figs
            saveas(abs_fig, fullfile(outdir, 'figures', ...
                strcat('pc_ode_', run_mode, '_abs.png')))
        end

        % initialize plot - relative amounts
        rel_fig = figure("Name","Pulse-chase time course plot - relative amounts");
        
        % Extracellular sFlt1 (X) relative amount
        subplot(2,1,1);
        plot(pc_sim_time, pc_sim_X_frac, 'LineWidth', 3)
        hold on
        scatter(jung_v2.time, jung_v2.media, 'filled')
        hold off
        xline(-pulse_length, '-k', "Start pulse", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Start chase", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Extracellular S Relative Amount (X)")
        xlabel("Time (h)")
        ylabel("X (normalized to 8h)")
        xlim([-1, 10]);
        ylim([-.05, 1.05]);
        
        % Intracellular sFlt1 (Si) amount
        subplot(2,1,2)
        plot(pc_sim_time, pc_sim_I_fc, 'LineWidth', 3)
        hold on
        scatter(jung_v2.time, jung_v2.lysate, 'filled')
        hold off
        xline(-pulse_length, '-k', "Start pulse", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Start chase", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Intracellular S Relative Amount (I)")
        xlabel("Time (h)")
        ylabel("I (normalized to 0h)")
        xlim([-1, 10]);

        % save figure - relative amounts of species
        if save_figs
            saveas(rel_fig, fullfile(outdir, 'figures', ...
                strcat('pc_ode_', run_mode, '_rel.png')))
        end

    case 'extended_pulse_chase'

        num_param_sets = 1;

        % lengthen pulse and chase to see long-term dynamics
        pulse_length = 24;   % pulse = 24 hours
        chase_length = 10*24;      % chase = 10 days
        time_interval = 10/60;   % record results every 10 minutes
        pulse_tspan = -pulse_length:time_interval:0;    % set pulse time span vector
        chase_tspan = 0:time_interval:chase_length;     % set chase time span vector

        module_hjk_opt_setup;

        % save simulation input parameters before running
        if save_params
            save(fullfile(outdir, 'inputs', 'pulse_chase_input.mat'), 'num_cells', 'pulse_length', ...
                'chase_length', 'time_interval', 'sp', 'p', 'conv_factor_ngml');
        end

        [pc_sim_time, pc_sim_Y, pc_sim_X_frac, pc_sim_I_fc] = sim_pulse_chase_ode(sp, ...
            p, pulse_length, chase_length, time_interval, conv_factor_ngml);
    
        % save time course, parameter values, options
        if save_data
            save(fullfile(outdir, 'results', 'pc_sim.mat'), ...
                'pc_sim_time', 'pc_sim_Y', 'pc_sim_X_frac', 'pc_sim_I_fc', ...
                'pulse_length', 'chase_length', 'time_interval');
        end

        % initialize plot - absolute amounts
        abs_fig = figure("Name","Pulse-chase time course plot - absolute amounts");

        % Extracellular sFlt1 (X) absolute amount
        subplot(2,1,1);
        plot(pc_sim_time, pc_sim_Y(:,sp.X), 'LineWidth', 3)
        xlim([-1.5*pulse_length, chase_length])
        xline(-pulse_length, '-k', "Start pulse", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Start chase", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Extracellular sFLT1 Amount (X)")
        xlabel("Time (h)")
        ylabel("Amount of X (ng/mL)")
        
        % Intracellular sFlt1 (I) absolute amount
        subplot(2,1,2)
        plot(pc_sim_time, pc_sim_Y(:,sp.I), 'LineWidth', 3)
        xlim([-1.5*pulse_length, chase_length])
        xline(-pulse_length, '-k', "Start pulse", 'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5)
        xline(0, '-k', "Start chase", 'LabelHorizontalAlignment', 'right', 'LineWidth', 1.5)
        title("Intracellular sFLT1 Amount (I)")
        xlabel("Time (h)")
        ylabel("Amount of I (#/cell)")


    case 'opt_1'

        rng(081423);

        tic

        %% SET UP FOR OPTIMIZATION

        % set cost options, number of parameter sets, initial parameter values;
        % set optimization options and bounds; initialize output arrays
        module_hjk_opt_setup;

        i = 1;

        % extract initial parameter values from param_mat
        %k_init = [p.alpha, p.beta, p.gamma, p.delta];  % debug
        k_init = param_mat(1,:);
        for j = 1:length(param_names)
            p.(param_names{j}) = k_init(j);
        end

        % define anonymous cost function to pass additional arguments
        cost_fxn = @(k_init) hjk_combo_cost(k_init, kinghorn.time, kinghorn.X_frac_24h, ...            
            kinghorn.I_fc_0h, jung_v2.time, jung_v2.media, jung_v2.lysate, hornig.basal_time, ...
            hornig.basal_Sx_ng_ml, sp, p, param_names, conv_factor_ngml,  pulse_length, ...
            chase_length, secr_end_time, time_interval, jung_cost_option, ...
            hornig_cost_option, kinghorn_cost_option);

        %% NONLINEAR LEAST SQUARES OPTIMIZATION

        % run nonlinear least squares optimization
        [optimal(i,:), cost(i), ~, exit_flag(i)] = ...
            lsqnonlin(cost_fxn, k_init, lb, ub, lsq_options);

        % set parameters to optimized values
        for j = 1:length(param_names)
            p.(param_names{j}) = optimal(i,j);
        end

        %% RUN PULSE-CHASE CASE WITH OPTIMAL PARAMETER SET
        module_jung_run_optimal;    

        module_constitutive_run_optimal;
           
        %% CALCULATE COST COMPONENTS

        jung_I_cost_ssd(i) = sum(jung_I_cost_detail(i,:).^2, "omitnan");
        jung_X_cost_ssd(i) = sum(jung_X_cost_detail(i,:).^2, "omitnan");    
        hornig_X_cost_ssd(i) = sum(hornig_X_cost_detail(i,:).^2, "omitnan");
        kinghorn_X_cost_ssd(i) = ...
            sum(kinghorn_X_cost_detail(i,:).^2, "omitnan");
        kinghorn_I_cost_ssd(i) = ...
            sum(kinghorn_I_cost_detail(i,:).^2, "omitnan");

        % test that cost from lsqnonlin matches calculated cost of fit
        % (error of less than cost/1000)
        cost_error = cost(i) - jung_I_cost_ssd(i) - ...
                jung_X_cost_ssd(i) - hornig_X_cost_ssd(i) - ...
                kinghorn_X_cost_ssd(i) - kinghorn_I_cost_ssd(i);
        if cost_error / cost(i) > .001
            warning("Cost from lsqnonlin does not match cost observed with " + ...
                "optimization run %d. \n", i);

            cost_mismatch(i) = 1;
        end


        % save optimal parameters and time course data for this cost/decay option set
        if save_data
            save(fullfile(outdir, 'results', 'hjk_opt_1.mat'), ...
                'optimal', 'cost', 'exit_flag', 'infinite_sadness', 'cost_mismatch', ...
                'jung_T_all', 'jung_Y_all', 'jung_X_frac_all', 'jung_I_fc_all', ...
                'jung_X_cost_detail', 'jung_X_cost_ssd', 'jung_I_cost_detail', ...
                'jung_I_cost_ssd', 'hornig_T_all', 'hornig_Y_all', ...
                'hornig_X_cost_detail', 'hornig_X_cost_ssd', ...
                'kinghorn_X_cost_ssd', 'kinghorn_X_cost_detail', ...
                'kinghorn_I_cost_ssd', 'kinghorn_I_cost_detail');

        end

        toc



    case 'opt_100'

        rng(081423);

        tic

        % set cost options, number of parameter sets, initial parameter values;
        % set optimization options and bounds; initialize output arrays
        module_hjk_opt_setup;
       

        % optimize model starting from n=100 random parameter sets
        for i = 1:num_param_sets
          
            % update counter
            fprintf('Optimizing model %d of %d...\n', i, num_param_sets);
    
            % set parameters to the ith set
            k_init = param_mat(i,:);
            for j = 1:length(param_names)
                p.(param_names{j}) = k_init(j);
            end

            % define anonymous cost function to pass additional arguments
            cost_fxn = @(k_init) hjk_combo_cost(k_init, kinghorn.time, kinghorn.X_frac_24h, ...            
                kinghorn.I_fc_0h, jung_v2.time, jung_v2.media, jung_v2.lysate, ...
                hornig.basal_time, hornig.basal_Sx_ng_ml, sp, p, param_names, ...
                conv_factor_ngml,  pulse_length, chase_length, secr_end_time, ...
                time_interval, jung_cost_option, hornig_cost_option, kinghorn_cost_option);

            % run nonlinear least squares optimization
            try
                [optimal(i,:), cost(i), ~, exit_flag(i)] = ...
                    lsqnonlin(cost_fxn, k_init, lb, ub, lsq_options);
            catch
                % reports when parameter sets fail to converge
                infinite_sadness(i) = 1;
                fprintf("Failed parameter sets: %d \n", sum(infinite_sadness));
                continue
            end

            % set parameters to optimized values
            for j = 1:length(param_names)
                p.(param_names{j}) = optimal(i,j);
            end

            try
                % run pulse-chase secretion case with optimal parameters
                module_jung_run_optimal;

                % run constitutive secretion case with optimal parameters
                module_constitutive_run_optimal;

            catch
                continue;
            end

            % calculate total cost of fit to each dataset
            jung_I_cost_ssd(i) = sum(jung_I_cost_detail(i,:).^2, "omitnan");
            jung_X_cost_ssd(i) = sum(jung_X_cost_detail(i,:).^2, "omitnan");    
            hornig_X_cost_ssd(i) = sum(hornig_X_cost_detail(i,:).^2, "omitnan");
            kinghorn_X_cost_ssd(i) = ...
                sum(kinghorn_X_cost_detail(i,:).^2, "omitnan");
            kinghorn_I_cost_ssd(i) = ...
                sum(kinghorn_I_cost_detail(i,:).^2, "omitnan");

            % test that cost from lsqnonlin matches calculated cost of fit
            % (error of less than cost/1000)
            cost_error = cost(i) - jung_I_cost_ssd(i) - ...
                    jung_X_cost_ssd(i) - hornig_X_cost_ssd(i) - ...
                    kinghorn_X_cost_ssd(i) - kinghorn_I_cost_ssd(i);
            if cost_error / cost(i) > .001
                warning("Cost from lsqnonlin does not match cost observed with " + ...
                    "optimization run %d. \n", i);

                cost_mismatch(i) = 1;
            end


        end



        disp("Saving optimized parameters and time courses, if desired...")


        % save optimal parameters and time course data for this cost/decay option set
if save_data
    save(fullfile(outdir, 'results', 'hjk_ode_opt.mat'), ...
        'optimal', 'cost', 'exit_flag', 'infinite_sadness', 'cost_mismatch', ...
        'jung_T_all', 'jung_Y_all', 'jung_X_frac_all', 'jung_I_fc_all', ...
        'jung_X_cost_detail', 'jung_X_cost_ssd', 'jung_I_cost_detail', ...
        'jung_I_cost_ssd', 'hornig_T_all', 'hornig_Y_all', ...
        'hornig_X_cost_detail', 'hornig_X_cost_ssd', ...
        'kinghorn_X_cost_ssd', 'kinghorn_X_cost_detail', ...
        'kinghorn_I_cost_ssd', 'kinghorn_I_cost_detail');

end
        
toc
        
end