%% Sensitivity analysis of sFLT1 DDE model optimized to Kinghorn+Jung+Hornig data

clear;
close all;

% set run mode
% base_const = constitutive secretion with base parameters
% local_sens_const = local (x%) sensitivity analysis of constitutive case parameters 
% sens_vary_inh_genchem = univariate sensitivity analysis - 
%       percent inhibition of individual parameters
% sens_logscale_sqrt10 = univariate sensitivity analysis -
%       scale individual parameters by powers of sqrt(10)
% explore_c1 = Vary alpha and beta while keeping c1 = alpha*beta constant
% explore_c2 = Vary beta and gamma while keeping c2 = beta+gamma constant
% explore_c1c2 = Vary alpha, beta, and gamma while keeping c1 = alpha*beta and c2 = beta+gamma constant
% sens_vary_base_inh_genchem = 
run_mode = "sens_vary_base_inh_genchem";

% Announce which case will run
disp("Running case " + run_mode + "...")

% Data saving options
save_params = 0;    % Option to save input parameters; 1 = yes, 0 = no
save_data = 0;  % Option to save data; 1 = yes, 0 = no
save_figs = 0;  % Option to save figures; 1 = yes, 0 = no

% set species, parameters, seed, output directory
module_setup;

switch run_mode

    case "base_const"

        % simulate constitutive results with optimal parameter set to store time course output
        [constitutive_sim_time, constitutive_sim_Y] = sim_secr_dde(sp, p, ...
            "on", mat_delay, internalize_case, secr_end_time, time_interval, ...
            conv_factor_ngml, status_message);

        % Extracellular sFlt1 (X) amount, ng/mL
        fig = figure;
        subplot(2,1,1);
        plot(constitutive_sim_time, constitutive_sim_Y(sp.X,:), 'LineWidth', 3)
        hold on;
        scatter(hornig.basal_time, hornig.basal_Sx_ng_ml)
        hold off;
        title("Extracellular Amount (X)")
        xlabel("Time (h)")
        ylabel("Amount (ng/ml)")
        
        % Intracellular sFlt1 (I) amount
        subplot(2,1,2)
        plot(constitutive_sim_time, constitutive_sim_Y(sp.I, :), 'LineWidth', 3)
        title("Intracellular Amount (I)")
        xlabel("Time (h)")
        ylabel("Amount (#/cell)")

        % save figure - absolute amounts of species
        if save_figs
            saveas(fig, fullfile(outdir, 'figures', ...
                strcat(run_mode, '.png')))
        end

        if save_data
            save base_const constitutive_sim_time constitutive_sim_Y param_names sp p
        end


    case "local_sens_const"

        % specify fractional change for parameters (ex .1 = +10%)
        frac_change = .1;

        % store initial p so it can be reset for new iterations
        p_init = p;

         % no chase - production constitutively on
        prod_case = "on";
        
        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p', 'p_init', 'param_names', 'conv_factor_ngml', 'frac_change');
        end
        
        % run simulation to find base case:
        % (run to steady state followed by media change and 72h followup)
        [sim_time, sim_y] = sim_secr_dde(sp, p, "on", mat_delay, internalize_case, ...
            secr_end_time, time_interval, conv_factor_ngml, status_message);

        % find baseline 72h values
        out = calc_baseline_outs(sim_time, sim_y, sp, p, conv_factor_ngml);

        % initialize arrays to store 72h output values and sensitivities after adjustments
        out = init_out_sens_arrays(out, length(param_names), 1);

        % test sensitivity of system to each parameter
        for i = 1:length(param_names)
         
            % reset parameters to initial values
            p = p_init;

            % adjust parameter i by desired % change
            p.(param_names{i}) = p.(param_names{i}) .* (1+frac_change);

            % run simulation to steady state followed by media change and 72h followup
            [sim_time, sim_y] = sim_secr_dde(sp, p, "on", mat_delay, internalize_case, ...
                secr_end_time, time_interval, conv_factor_ngml, status_message);
            

            % calculate and store output values and sensitivities for run i
            out = store_out_val_raw_set(out, i, 1, sim_time, sim_y, sp, p, conv_factor_ngml);
            out = store_out_val_rel_sens_set(out, i, 1, frac_change);

        end

        % save outputs
        if save_data
            save(fullfile(outdir, 'results', 'local_sens_const.mat'), 'out', 'param_names')
        end

    case "sens_vary_inh_genchem"

        % select range of inhibition %s
        inh_frac = linspace(.01, 1, 100);

        % specify parameter modification type
        % genetic: apply modification before running system to steady state
        % chemical: run to steady state, then apply modification at t=0
        genetic_chemical = "chemical";

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);

        % store original parameter values 
        p_init = p;

        % introduce p_inh: for now, it's the same as p;
        % will modify in later runs
        p_inh = p;

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'conv_factor_ngml');
        end

        %% Find base case and initialize output arrays

        % run simulation to find base case (recall here p_inh == p)
        [sim_time, base_sim_y] = sim_secr_dde_inh_t0(sp, p, p_inh, secr_end_time,...
                time_interval, conv_factor_ngml, status_message);

        % extract X and I time courses
        X_base_tc = base_sim_y(sp.X,:);
        I_base_tc = base_sim_y(sp.I,:);

        % find baseline 72h values
        out = calc_baseline_outs(sim_time, base_sim_y, sp, p, conv_factor_ngml);

        % initialize arrays to store 72h output values and sensitivities after adjustments
        % dimensions: (parameter varied, percent inhibition)
        out = init_out_sens_arrays(out, length(param_names), length(inh_frac));

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % initialize container for storing input parameter values
        % matrix p_inh_all(:,n) gives parameter values for nth parameter set
        p_inh_all = ones(length(param_names), length(param_names) .* length(inh_frac)) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        %T_all = ones(num_sim_pts, length(param_names) .* length(inh_case)) * NaN; 
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), length(param_names) .* length(inh_frac)) * NaN;

        %% Run experimental cases
        % find fold change in X_timept, I_timept when decreasing param_names{i} by inh_frac (j)
        for i = 1:length(param_names)

            for j = 1:length(inh_frac)

                % reset parameters to initial values
                p = p_init;
                p_inh = p_init;

                % adjust parameter i by desired % inhibition
                p_inh.(param_names{i}) = p_init.(param_names{i}) .* (1-inh_frac(j));

                % calculate index of this parameter set
                out_index = (i-1) .* length(inh_frac) + j;

                % store p_inh values for data analysis
                for h = 1:length(param_names)
                    p_inh_all(h, out_index) = p_inh.(param_names{h});
                end

                %% branch point - genetic inh or chemical inh?
                switch genetic_chemical

                    % genetic: perturbation applied at beginning before
                    % reaching steady state
                    case "genetic"

                        % run both parts of simulation with p_inh
                        try
                            [~, sim_y] = sim_secr_dde_inh_t0(sp, p_inh, p_inh, secr_end_time,...
                                time_interval, conv_factor_ngml, status_message);
                        catch
                            fprintf("Failed parameter sets: %s x %d \n", param_names{i}, 1 - inh_frac(j))
                            continue
                        end
        
                    % chemical: perturbation applied at t=0
                    case "chemical"

                        % run first part with p, then second with p_inh
                        try
                            [~, sim_y] = sim_secr_dde_inh_t0(sp, p, p_inh, secr_end_time,...
                                time_interval, conv_factor_ngml, status_message);
                        catch
                            fprintf("Failed parameter sets: %s x %d \n", param_names{i}, 1 - inh_frac(j))
                            continue
                        end
        
                end

                % store phenomenological parameters (72h values, AUCs, etc)
                out = store_out_val_raw_set(out, i, j, sim_time, sim_y, sp, p, conv_factor_ngml);
                out = store_out_val_rel_sens_set(out, i, j, inh_frac(j));

                % store time courses in a single data structure
                Y_all(:,:,out_index) = sim_y';

            end

        end

        % save outputs
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('sens_vary_inh_genchem_', genetic_chemical, '.mat')), ... 
            'out', 'param_names', 'inh_frac', 'genetic_chemical', 'Y_all', 'sim_time', 'base_sim_y', ...
            'X_base_tc', "I_base_tc")
        end

case "sens_logscale_sqrt10"        

        % specify parameter modification type
        % genetic: apply modification before running system to steady state
        % chemical: run to steady state, then apply modification at t=0
        genetic_chemical = "genetic";

        % set scaling factors to try
        factors_to_try = logspace(-2, 2, 9);  % sqrt(10) scaling from .01 to 100

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);

        % store original parameter values 
        p_init = p;

        % introduce p_mod: for now, it's the same as p;
        % will modify in later runs
        p_mod = p;

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'conv_factor_ngml');
        end

        %% Find base case and initialize output arrays

        % run simulation to find base case (recall here p_mod == p)
        [sim_time, base_sim_y] = sim_secr_dde_inh_t0(sp, p, p_mod, secr_end_time,...
                time_interval, conv_factor_ngml, status_message);

        % extract X and I time courses
        X_base_tc = base_sim_y(sp.X,:);
        I_base_tc = base_sim_y(sp.I,:);

        % find baseline 72h values
        out = calc_baseline_outs(sim_time, base_sim_y, sp, p, conv_factor_ngml);

        % initialize arrays to store 72h output values and sensitivities after adjustments
        % dimensions: (parameter varied, percent inhibition)
        out = init_out_raw_arrays(out, length(param_names), length(factors_to_try));
        out = init_out_sens_arrays(out, length(param_names), length(factors_to_try));

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % initialize container for storing input parameter values
        % matrix p_mod_all(:,n) gives parameter values for nth parameter set
        p_mod_all = ones(length(param_names), length(param_names) .* length(factors_to_try)) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        %T_all = ones(num_sim_pts, length(param_names) .* length(inh_case)) * NaN; 
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), length(param_names) .* length(factors_to_try)) * NaN;

        %% Run experimental cases
        % find fold change in X_timept, I_timept when multiplying param_names{i} by factors_to_try(j)
        for i = 1:length(param_names)

            for j = 1:length(factors_to_try)

                % reset parameters to initial values
                p = p_init;
                p_mod = p_init;

                % adjust parameter i by desired % inhibition
                p_mod.(param_names{i}) = p_init.(param_names{i}) .* factors_to_try(j);

                % calculate index of this parameter set
                out_index = (i-1) .* length(factors_to_try) + j;

                % store p_mod values for data analysis
                for h = 1:length(param_names)
                    p_mod_all(h, out_index) = p_mod.(param_names{h});
                end

                %% branch point - genetic or chemical perturbation?
                switch genetic_chemical

                    % genetic: perturbation applied at beginning before
                    % reaching steady state
                    case "genetic"

                        % run both parts of simulation with p_mod
                        try
                            [~, sim_y] = sim_secr_dde_inh_t0(sp, p_mod, p_mod, secr_end_time,...
                                time_interval, conv_factor_ngml, status_message);
                        catch
                            fprintf("Failed parameter sets: %s x %d \n", param_names{i}, factors_to_try(j))
                            continue
                        end
        
                    % chemical: perturbation applied at t=0
                    case "chemical"

                        % run first part with p, then second with p_mod
                        try
                            [~, sim_y] = sim_secr_dde_inh_t0(sp, p, p_mod, secr_end_time,...
                                time_interval, conv_factor_ngml, status_message);
                        catch
                            fprintf("Failed parameter sets: %s x %d \n", param_names{i}, factors_to_try(j))
                            continue
                        end
        
                end

                % store phenomenological parameters (72h values, AUCs, etc)
                out = store_out_val_raw_set(out, i, j, sim_time, sim_y, sp, p, conv_factor_ngml);
                
                % store relative sensitivities vs baseline case
                out = store_out_val_rel_sens_set(out, i, j, factors_to_try(j));

                % store time courses in a single data structure
                Y_all(:,:,out_index) = sim_y';

            end

        end

        % save outputs
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('sens_logscale_sqrt10_', genetic_chemical, '.mat')), ... 
            'out', 'param_names', 'factors_to_try', 'p_mod_all', 'genetic_chemical', 'Y_all', ...
            'sim_time', 'base_sim_y', 'X_base_tc', "I_base_tc")
        end

case "explore_c1"
        % definition: c1 = alpha * beta
        p_init = p;  % store original parameter values (medians from optimization)
        c1 = p_init.alpha * p_init.beta;
        %c2 = p_init.beta + p_init.gamma;    % used to calculate range for beta

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);


        % designate range of factors for multiplying (powers of 2 from 1/8 to 8)
        factors = 2 .^ (-3:3);

        % no chase - production constitutively on
        prod_case = "on";

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'conv_factor_ngml', 'factors');
        end

        % initialize containers for phenomenological parameters
        % dimensions: (beta/gamma combinations, 1)
        out = init_out_raw_arrays(struct(), length(factors), 1);

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % initialize container for storing input parameter values
        % matrix p_mod_all(:,n) gives parameter values for nth parameter set
        p_mod_all = ones(length(param_names), length(factors)) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), length(factors)) * NaN;

        % test different distributions of c1 across (alpha, beta)
        for i = 1:length(factors)
            p_mod = p_init;
            p_mod.alpha = p_init.alpha .* factors(i);
            p_mod.beta = p_init.beta ./ factors(i);
            
            % store p_mod values for data analysis
            for h = 1:length(param_names)
                p_mod_all(h, i) = p_mod.(param_names{h});
            end
            
            % run simulation to steady state followed by media change and 72h followup
            [sim_time, sim_y] = sim_secr_dde(sp, p_mod, prod_case, mat_delay, internalize_case, ...
                secr_end_time, time_interval, conv_factor_ngml, status_message);

            % calculate and store output values for run i
            out = store_out_val_raw_set(out, i, 1, sim_time, sim_y, sp, p_mod, conv_factor_ngml);

            % store time courses in a single data structure
            Y_all(:,:,i) = sim_y';

        end

        % save output
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('explore_c1.mat')), ...
                'out', 'param_names', 'p_mod_all', 'sim_time',  ...
                'Y_all')
        end


    case "explore_c2"
        % definition: c2 = beta + gamma
        p_init = p;  % store original parameter values (medians from optimization)

        c2 = p_init.beta + p_init.gamma;

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);

        % store original parameter values 
        p_init = p;

        % designate number of simulations to test
        num_parts = 6;  % number of sims will be num_parts + 1
        chunk_size = c2/num_parts;

        % no chase - production constitutively on
        prod_case = "on";

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'conv_factor_ngml', 'num_parts', 'chunk_size');
        end

        % initialize containers for phenomenological parameters
        % dimensions: (beta/gamma combinations, 1)
        out = init_out_raw_arrays(struct(), num_parts + 1, 1);

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % initialize container for storing input parameter values
        % matrix p_mod_all(:,n) gives parameter values for nth parameter set
        p_mod_all = ones(length(param_names), num_parts + 1) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), num_parts + 1) * NaN;

        % test different distributions of c2 across (beta, gamma)
        % example: if num_parts = 3,
        % test (b,g) sets = (0, c2), (c2/3, 2c2/3), (2c2/3, c2/3), (c2, 0)
        for i = 1:(num_parts + 1)
            p_mod = p;
            p_mod.beta = (i-1) * chunk_size;
            p_mod.gamma = c2 - p_mod.beta;
            
            % store p_mod values for data analysis
            for h = 1:length(param_names)
                p_mod_all(h, i) = p_mod.(param_names{h});
            end

            % run simulation to steady state followed by media change and 72h followup
            [sim_time, sim_y] = sim_secr_dde(sp, p_mod, prod_case, mat_delay, internalize_case, ...
                secr_end_time, time_interval, conv_factor_ngml, status_message);

            % calculate and store output values for run i
            out = store_out_val_raw_set(out, i, 1, sim_time, sim_y, sp, p_mod, conv_factor_ngml);

            % store time courses in a single data structure
            Y_all(:,:,i) = sim_y';

        end

        % save output
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('explore_c2_', int2str(num_parts), '_parts.mat')), ...
                'out', 'param_names', 'p_mod_all', 'sim_time',  ...
                'Y_all')
        end

    case "explore_c1c2"

        % definitions: c1 = alpha * beta, c2 = beta + gamma
        p_init = p;  % store original parameter values (medians from optimization)

        c1 = p_init.alpha * p_init.beta;
        c2 = p_init.beta + p_init.gamma;

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);

        % store original parameter values 
        p_init = p;

        % generate num_steps values for beta from 0 to c2, excluding 0
        num_steps = 100;
        betas = linspace(0, c2, num_steps + 1);
        betas = betas(2:end);

        % no chase - production constitutively on
        prod_case = "on";

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'conv_factor_ngml', 'num_steps');
        end

        % initialize containers for phenomenological parameters
        % dimensions: (alpha,beta,gamma) sets, 1)
        out = init_out_raw_arrays(struct(), length(betas), 1);

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % initialize container for storing input parameter values
        % matrix p_mod_all(:,n) gives parameter values for nth parameter set
        p_mod_all = ones(length(param_names), length(betas)) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), length(betas)) * NaN;

        % test different distributions of c2 across (beta, gamma)
        % example: if num_parts = 3,
        % test (b,g) sets = (0, c2), (c2/3, 2c2/3), (2c2/3, c2/3), (c2, 0)
        for i = 1:length(betas)
            p_mod = p;
            p_mod.beta = betas(i);
            p_mod.alpha = c1 / p_mod.beta;
            p_mod.gamma = c2 - p_mod.beta;
            
            % store p_mod values for data analysis
            for h = 1:length(param_names)
                p_mod_all(h, i) = p_mod.(param_names{h});
            end
            
            % run simulation to steady state followed by media change and 72h followup
            [sim_time, sim_y] = sim_secr_dde(sp, p_mod, prod_case, mat_delay, internalize_case, ...
                secr_end_time, time_interval, conv_factor_ngml, status_message);

            % calculate and store output values for run i
            out = store_out_val_raw_set(out, i, 1, sim_time, sim_y, sp, p_mod, conv_factor_ngml);

            % store time courses in a single data structure
            Y_all(:,:,i) = sim_y';

        end

        % save output
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('explore_c1c2_', int2str(num_steps), '_steps.mat')), ...
                'out', 'param_names', 'p_mod_all', 'sim_time',  ...
                'Y_all')
        end


    case "sens_vary_base_inh_genchem"

        p_init = p;  % store original parameter values (medians from optimization)

        % definitions: c1 = alpha * beta, c2 = beta + gamma
        c1 = p_init.alpha * p_init.beta;
        c2 = p_init.beta + p_init.gamma;

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);

        % generate num_steps values for beta from 0 to c2, excluding 0
        num_steps = 5;
        betas = linspace(0, c2, num_steps + 1);
        betas = betas(2:end);
        
        % initialize container for storing base parameter values
        % matrix p_base_all(:,h) gives the hth set of base parameter values
        p_base_all = ones(length(param_names), length(betas)) * NaN;
        
        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;
        num_sim_pts = length(opt_tspan);  % number of simulated time points

        % initialize container for base time courses and output values
        out = init_out_base_arrays(struct(), length(betas), num_sim_pts);

        % select range of inhibition %s
        inh_frac = linspace(.01, .99, 99);

        % calculate total number of runs to perform
        total_runs = length(betas) * length(param_names) * length(inh_frac);

        % initialize container for storing input parameter values
        % matrix p_inh_all(:,n) gives parameter values for nth parameter set
        p_inh_all = ones(length(param_names), total_runs) * NaN;

        % initialize containers for phenomenological parameters
        out = init_out_raw_arrays(out, total_runs, 1);

        % initialize arrays to store 72h output values and sensitivities after adjustments
            % dimensions: (parameter varied, percent inhibition)
        out = init_out_sens_arrays(out, total_runs, 1);
                
        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), total_runs) * NaN;

        % initialize p_inh: for now, it's the same as p;
        % will modify in later runs
        p_inh = p;

        % specify parameter modification type
        % genetic: apply modification before running system to steady state
        % chemical: run to steady state, then apply modification at t=0
        genetic_chemical = "chemical";

        % no chase - production constitutively on
        prod_case = "on";

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'betas', 'conv_factor_ngml');
        end

        % for each beta value, test inhibition of each parameter 
        % by each inh_frac
        for h = 1:length(betas)
            % set parameters
            p_base = p;
            p_base.beta = betas(h);
            p_base.alpha = c1 / p_base.beta;
            p_base.gamma = c2 - p_base.beta;
            
            % store p_base values for data analysis
            for k = 1:length(param_names)
                p_base_all(k, h) = p_base.(param_names{k});
            end
            
            % run baseline simulation to steady state followed by media change and 72h followup
            [sim_time, base_sim_y] = sim_secr_dde(sp, p_base, prod_case, mat_delay, internalize_case, ...
                secr_end_time, time_interval, conv_factor_ngml, status_message);

            % save baseline 72h output values
            out = calc_baseline_outs_diffbase(out, h, sim_time, ...
                base_sim_y, sp, p_base, conv_factor_ngml);


        %% Run experimental cases
        % find fold change in X_timept, I_timept when decreasing param_names{i} by inh_frac (j)
        for i = 1:length(param_names)

            for j = 1:length(inh_frac)

                % set parameters to current base values
                p = p_base;
                p_inh = p_base;

                % adjust parameter i by desired % inhibition
                p_inh.(param_names{i}) = p_base.(param_names{i}) .* (1-inh_frac(j));

                % calculate index of this parameter set
                out_index = (h-1)*length(param_names)*length(inh_frac) + ...
                    (i-1)*length(inh_frac) + j;

                % store p_inh values for data analysis
                for k = 1:length(param_names)
                    p_inh_all(k, out_index) = p_inh.(param_names{k});
                end

                %% branch point - genetic inh or chemical inh?
                switch genetic_chemical

                    % genetic: perturbation applied at beginning before
                    % reaching steady state
                    case "genetic"

                        % run both parts of simulation with p_inh
                        try
                            [~, sim_y] = sim_secr_dde_inh_t0(sp, p_inh, p_inh, secr_end_time,...
                                time_interval, conv_factor_ngml, status_message);
                        catch
                            fprintf("Failed parameter sets: %s x %d \n", param_names{i}, 1 - inh_frac(j))
                            continue
                        end
        
                    % chemical: perturbation applied at t=0
                    case "chemical"

                        % run first part with p, then second with p_inh
                        try
                            [~, sim_y] = sim_secr_dde_inh_t0(sp, p, p_inh, secr_end_time,...
                                time_interval, conv_factor_ngml, status_message);
                        catch
                            fprintf("Failed parameter sets: %s x %d \n", param_names{i}, 1 - inh_frac(j))
                            continue
                        end
        
                end

                % store phenomenological parameters (72h values, AUCs, etc)
                out = store_out_val_raw_set_diffbase(out, out_index, sim_time, sim_y, sp, p, conv_factor_ngml);
                out = store_out_val_rel_sens_set_diffbase(out, h, out_index, inh_frac(j));

                % store time courses in a single data structure
                Y_all(:,:,out_index) = sim_y';

                end
            end
        end

        % CHECK OUTPUTS

        % save simulation inputs 
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_base_all', 'betas', 'conv_factor_ngml');
        end

        % save outputs
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('sens_diffbase_vary_inh_genchem_', ...
                  genetic_chemical, '.mat')), ...
                'out', 'param_names', 'inh_frac', 'genetic_chemical',  ...
                'Y_all', 'sim_time')
        end


end




        
        

        

        


















%% CASES TO ADD (original code, may need updates)
% explore_bg = vary both beta and gamma from 0:1
switch run_mode

    case "ORIG_explore_bg"

        % specify parameter modification type
        % genetic: apply modification before running system to steady state
        % chemical: run to steady state, then apply modification at t=0
        genetic_chemical = "genetic";

        % specify parameter ranges to test

        bg_ranges = "by_pt1_to_7";

        % specify case for tau values
        % "jh_opt": median value from Jung-Hornig optimization
        % "1234": 1, 2, 3, 4
        % "by_2s": powers of 2 from 1/8 to 8
        % "pt5_to_1": seq(0.5,1, by = 0.5)
        % "pt75_to_pt8": seq(0.75,0.8, by = 0.001)
        tau_range = "by_2s";

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);
       
        % set to constitutive production
        prod_case = "on";

        % calculate a refined threshold with more significant digits
        % (logical)
        refine_threshold = 1;

        % specify range of beta and gamma values to try
        switch bg_ranges
            case "by_pt1"
                betas = 0:0.1:1;
                gammas = 0:0.1:1;
            case "by_pt01"
                betas = 0:0.01:1;
                gammas = 0:0.01:1;
            case "explore_pt8_pt9_region"
                betas = 0.4;
                gammas = 0.4:.001:0.5;
            case "by_pt1_to_7"
                betas = 0:0.1:7;
                gammas = 0:0.1:7;
        end

        switch tau_range
            case "jh_opt"
                taus = 1.76;
            case "hjk_opt"
                taus = 1.96;
            case "1234"
                taus = 1:4;
            case "by_2s"
                taus = 2.^(-3:3);
            case "pt5_to_1"
                taus = 0.5:.05:1;
            case "pt75_to_pt8"
                taus = 0.75:.001:0.8;
        end

        % store original parameter values 
        p_init = p;

        % introduce modified parameter values; for now, set to init values
        p_mod = p;

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'betas', 'gammas', 'taus', 'conv_factor_ngml');
        end

        % initialize arrays to store 72h output values and sensitivities after adjustments
        % dimensions: (beta length, gamma length, tau length)
        %out = struct();
        out.I_ss = ones(length(betas), length(gammas), length(taus)) * NaN;
        out.t_crash = ones(length(betas), length(gammas), length(taus)) * NaN;

        % initialize containers for beta+gamma thresholds at each tau value
        bg_threshold_ub_at_tau = ones(length(taus), 1) * NaN;
        bg_threshold_lb_at_tau = ones(length(taus), 1) * NaN;

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % matrix p_mod_all(:,n) gives parameter values for nth parameter
        % set, n = out_index
        p_mod_all = ones(length(param_names), ...
            length(taus) .* length(betas) * length(gammas)) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), ...
            length(taus) .* length(betas) .* length(gammas)) * NaN;
                        
        for h = 1:length(taus)
            % note that threshold has not yet been found (or bounded) for this tau value
            bg_threshold_ub = NaN;
            bg_threshold_lb = NaN;


            for i = 1:length(betas)
                for j = 1:length(gammas)
                    % set current parameter values
                    p_mod.tau = taus(h);
                    p_mod.beta = betas(i);
                    p_mod.gamma = gammas(j);

                    % calculate index of this parameter set
                    out_index = (h-1) * length(betas) * length(gammas) + ...
                        (i-1) * length(betas) + j;

                    % store p_mod values for data analysis
                    for x = 1:length(param_names)
                        p_mod_all(x, out_index) = p_mod.(param_names{x});
                    end
                
                    % Simulation to steady state
                    try
                        [sim_time, sim_y] = sim_secr_dde_to_ss(sp, p_mod, prod_case, secr_end_time,...
                            time_interval, conv_factor_ngml, status_message);
    
                        % store phenomenological parameters (72h values, AUCs, etc)
                        % out = store_out_val_raw_set(out, i, j, sim_time, sim_y, sp, p, conv_factor_ngml);
                    
                        % store values of I and X at 72h
                        %out.X_ss(i,j) = sim_y(sp.X, end);    % ng/mL
                        out.I_ss(i,j,h) = sim_y(sp.I, end);    % #/cell

                        % find t_crash?
                        out.t_crash(i,j,h) = sim_time(dsearchn(sim_y(sp.I,:)', -1));

                        % store Y_all? does not yet work
%                        Y_all(1:size(sim_y,2),:,out_index) = sim_y';
    
                        % if the simulation failed bc I < 0, find b+g
                        % and update failure threshold if needed 
                        if sim_y(sp.I, end) < 0
    
                            % calculate beta + gamma
                            bg = betas(i) + gammas(j);
    
                            % initialize b+g threshold if needed
                            % (note calculated once for each tau)
                            if isnan(bg_threshold_ub)
                                bg_threshold_ub = bg;
                            end
                            
                            % update b+g threshold if needed
                            if bg < bg_threshold_ub
                                bg_threshold_ub = bg;
                            end
                        end
    
                    % note any failed parameter sets
                    catch
                        fprintf("Failed parameter sets: %d x %d x %d\n", ...
                            betas(i), gammas(j), taus(h))
                        continue
                    end
                end
            end
                    
            % refine beta+gamma threshold for this tau value 
            % to 2 additional significant digits    
            if refine_threshold && ~isnan(bg_threshold_ub)
                new_range_size = betas(2) - betas(1);
                bg_threshold_lb = bg_threshold_ub - new_range_size;
                new_step_size = new_range_size / 100;
        
                % since all b,g w same b+g have same results,
                % fix 1 and vary other: b+g = lb, b+g+range = ub
                % set beta to lb/2; vary gamma from 0.5*lb to lb/2+range
                new_beta = bg_threshold_lb / 2;
                new_gammas = new_beta:new_step_size:(new_beta + new_range_size);
        
                % Run experimental cases for grid of beta and gamma
                [bg_threshold_ub, ~, ~] = find_bg_threshold(new_beta,  ...
                    new_gammas, bg_threshold_ub, sp, p_mod, prod_case, secr_end_time, ...
                    time_interval, conv_factor_ngml, status_message,  param_names, p_mod_all);
        
                bg_threshold_lb = bg_threshold_ub - new_step_size;

            end

            % store beta+gamma threshold bounds for this value of tau
            bg_threshold_ub_at_tau(h) = bg_threshold_ub;
            bg_threshold_lb_at_tau(h) = bg_threshold_lb;
                
        end

        % save output
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('explore_bg_', bg_ranges, '_tau_', tau_range, '.mat')), ...
                'out', 'param_names', 'p_mod_all', 'bg_threshold_ub_at_tau',  ...
                'bg_threshold_lb_at_tau')
        end
        
    case "ORIG_explore_bg_alpha"

        % specify parameter modification type
        % genetic: apply modification before running system to steady state
        % chemical: run to steady state, then apply modification at t=0
        genetic_chemical = "genetic";

        % specify parameter ranges to test

        bg_ranges = "by_pt1";

        % specify case for alpha values
        % "jh_opt": median value from Jung-Hornig optimization
        % "by_2s": scale base value by powers of 2 from 1/8 to 8
        % "by_sqrt10s: scale base value by powers of sqrt(10) from .01x too 
        alpha_range = "by_sqrt10s";

        % set all initial species values to 0
        y0 = zeros(length(fieldnames(sp)), 1);
       
        % set to constitutive production
        prod_case = "on";

        % calculate a refined threshold with more significant digits
        % (logical)
        refine_threshold = 1;

        % specify range of beta and gamma values to try
        switch bg_ranges
            case "by_pt1"
                betas = 0:0.1:1;
                gammas = 0:0.1:1;
            case "by_pt01"
                betas = 0:0.01:1;
                gammas = 0:0.01:1;
            case "explore_pt8_pt9_region"
                betas = 0.4;
                gammas = 0.4:.001:0.5;
            case "by_pt1_to_7"
                betas = 0:0.1:7;
                gammas = 0:0.1:7;
        end

        % store original parameter values 
        p_init = p;

        switch alpha_range
            case "jh_opt"
                alphas = p_init.alpha;
            case "by_2s"
                alphas = p_init.alpha .* 2 .^ (-3:3);
            case "by_sqrt10s"
                alphas =  p_init.alpha .* sqrt(10) .^ (-4:4);
        end
        
        % introduce modified parameter values; for now, set to init values
        p_mod = p;

        % save simulation input before running to steady state
        if save_params
            save(fullfile(outdir, 'inputs', 'input.mat'), 'num_cells', 'secr_end_time', ...
                'y0', 'sp', 'p_init', 'betas', 'gammas', 'alphas', 'conv_factor_ngml');
        end

        % initialize arrays to store 72h output values and sensitivities after adjustments
        % dimensions: (beta length, gamma length, tau length)
        %out = struct();
        out.I_ss = ones(length(betas), length(gammas), length(alphas)) * NaN;
        out.t_crash = ones(length(betas), length(gammas), length(alphas)) * NaN;

        % initialize containers for beta+gamma thresholds at each tau value
        bg_threshold_ub_at_alpha = ones(length(alphas), 1) * NaN;
        bg_threshold_lb_at_alpha = ones(length(alphas), 1) * NaN;

        % set time span for reporting output
        opt_tspan = 0:time_interval:secr_end_time;

        % matrix p_mod_all(:,n) gives parameter values for nth parameter
        % set, n = out_index
        p_mod_all = ones(length(param_names), ...
            length(alphas) .* length(betas) * length(gammas)) * NaN;

        % initialize output containers for simulation data
        % num_sim_pts: number of points in simulation output
        % vector T_all(:,n) gives time series for nth parameter set
        % matrix Y_all(:,:,n) gives simulation output matrix for nth param set
        num_sim_pts = length(opt_tspan);
        Y_all = ones(num_sim_pts, length(fieldnames(sp)), ...
            length(alphas) .* length(betas) .* length(gammas)) * NaN;
                        
        for h = 1:length(alphas)
            % note that threshold has not yet been found (or bounded) for this tau value
            bg_threshold_ub = NaN;
            bg_threshold_lb = NaN;


            for i = 1:length(betas)
                for j = 1:length(gammas)
                    % set current parameter values
                    p_mod.alpha = alphas(h);
                    p_mod.beta = betas(i);
                    p_mod.gamma = gammas(j);

                    % calculate index of this parameter set
                    out_index = (h-1) * length(betas) * length(gammas) + ...
                        (i-1) * length(betas) + j;

                    % store p_mod values for data analysis
                    for x = 1:length(param_names)
                        p_mod_all(x, out_index) = p_mod.(param_names{x});
                    end
                
                    % Simulation to steady state
                    try
                        [sim_time, sim_y] = sim_secr_dde_to_ss(sp, p_mod, prod_case, secr_end_time,...
                            time_interval, conv_factor_ngml, status_message);
    
                        % store phenomenological parameters (72h values, AUCs, etc)
                        % out = store_out_val_raw_set(out, i, j, sim_time, sim_y, sp, p, conv_factor_ngml);
                    
                        % store values of I and X at 72h
                        %out.X_ss(i,j) = sim_y(sp.X, end);    % ng/mL
                        out.I_ss(i,j,h) = sim_y(sp.I, end);    % #/cell

                        % find t_crash?
                        out.t_crash(i,j,h) = sim_time(dsearchn(sim_y(sp.I,:)', -1));

                        % store Y_all? does not yet work
%                        Y_all(1:size(sim_y,2),:,out_index) = sim_y';
    
                        % if the simulation failed bc I < 0, find b+g
                        % and update failure threshold if needed 
                        if sim_y(sp.I, end) < 0
    
                            % calculate beta + gamma
                            bg = betas(i) + gammas(j);
    
                            % initialize b+g threshold if needed
                            % (note calculated once for each tau)
                            if isnan(bg_threshold_ub)
                                bg_threshold_ub = bg;
                            end
                            
                            % update b+g threshold if needed
                            if bg < bg_threshold_ub
                                bg_threshold_ub = bg;
                            end
                        end
    
                    % note any failed parameter sets
                    catch
                        fprintf("Failed parameter sets: %d x %d x %d\n", ...
                            betas(i), gammas(j), alphas(h))
                        continue
                    end
                end
            end
                    
            % refine beta+gamma threshold for this tau value 
            % to 2 additional significant digits    
            if refine_threshold && ~isnan(bg_threshold_ub)
                new_range_size = betas(2) - betas(1);
                bg_threshold_lb = bg_threshold_ub - new_range_size;
                new_step_size = new_range_size / 100;
        
                % since all b,g w same b+g have same results,
                % fix 1 and vary other: b+g = lb, b+g+range = ub
                % set beta to lb/2; vary gamma from 0.5*lb to lb/2+range
                new_beta = bg_threshold_lb / 2;
                new_gammas = new_beta:new_step_size:(new_beta + new_range_size);
        
                % Run experimental cases for grid of beta and gamma
                [bg_threshold_ub, ~, ~] = find_bg_threshold(new_beta,  ...
                    new_gammas, bg_threshold_ub, sp, p_mod, prod_case, secr_end_time, ...
                    time_interval, conv_factor_ngml, status_message,  param_names, p_mod_all);
        
                bg_threshold_lb = bg_threshold_ub - new_step_size;

            else
                bg_threshold_lb = bg_threshold_ub - (betas(2) - betas(1));
            end

            % store beta+gamma threshold bounds for this value of tau
            bg_threshold_ub_at_alpha(h) = bg_threshold_ub;
            bg_threshold_lb_at_alpha(h) = bg_threshold_lb;
                
        end

        % save output
        if save_data
            save(fullfile(outdir, 'results', ...
                strcat('explore_bg_', bg_ranges, '_alpha_', alpha_range, '.mat')), ...
                'out', 'param_names', 'p_mod_all', 'bg_threshold_ub_at_alpha',  ...
                'bg_threshold_lb_at_alpha')
        end


    case "ORIG_random_plot_code_1"

        t = linspace(0, 1, 61); % every min for 1h
        plot(t, exp(-9*t), "LineWidth",3);
        hold on;
        plot(t, exp(-16*t), "LineWidth", 3);
        hold off;
        legend({"\kappa = 9", "\kappa = 16"}, "FontSize", 20);
        xlabel("Time (h)");
        ylabel("Fraction of initial secretion", "FontSize", 20);
        ax = gca;
        ax.FontSize = 20;

    case "ORIG_random_plot_code_2"

        plot(inh_frac, out.X_72(p_ind.beta, :) ./ out.X_72_0, ...
            "LineWidth",3)
        xlabel("Fraction inhibition of \beta")
        ylabel("X at 72h (fraction of control)")
        ax = gca;
        ax.FontSize = 20;

end

disp("Run complete!")





%% FUNCTIONS

function out = calc_baseline_outs(sim_time, sim_y, sp, p, conv_factor_ngml)
    % find baseline 72h output values
    
    % species amounts at 72h
    out.X_72_0 = sim_y(sp.X, end); % ng/mL
    out.I_72_0 = sim_y(sp.I, end); % #/cell
    % species AUC at 72h
    out.AUC_X_72_0 = trapz(sim_time, sim_y(sp.X,:)); % (ng * h)/mL
    out.AUC_I_72_0 = trapz(sim_time, sim_y(sp.I,:)); % (# * h)/cell
    % species AUMC at 72h
    out.AUMC_X_72_0 = trapz(sim_time, sim_time .* sim_y(sp.X,:)); % (ng * h^2)/mL
    out.AUMC_I_72_0 = trapz(sim_time, sim_time .* sim_y(sp.I,:)); % (# * h^2)/cell
    % species MRT at 72h - not sure this works without AUC_inf, AUMC_inf
    %out.MRT_X_72_0 = out.AUMC_X_72_0 / out.AUC_X_72_0; % h
    %out.MRT_I_72_0 = out.AUMC_I_72_0 / out.AUC_I_72_0; % h
    % analytic steady state for I, X
    out.I_ss_calc_0 = p.alpha / (p.beta + p.gamma);
    out.X_ss_calc_0 = p.alpha * p.beta * conv_factor_ngml/ (p.delta * (p.beta + p.gamma));
    % T50 for X_ss (time when X is halfway to analytic steady state
    out.T50_X_0 = sim_time(dsearchn(sim_y(sp.X,:)', out.X_ss_calc_0/2));

end

function out = init_out_raw_arrays(out, dim1, dim2)
    % initialize arrays to store 72h values and sensitivities after adjustments

    % arrays for 72h values
    out.X_72 = ones(dim1,dim2) * NaN;
    out.I_72 = ones(dim1,dim2) * NaN;
    out.AUC_X_72 = ones(dim1,dim2) * NaN;
    out.AUC_I_72 = ones(dim1,dim2) * NaN;
    out.AUMC_X_72 = ones(dim1,dim2) * NaN;
    out.AUMC_I_72 = ones(dim1,dim2) * NaN;
    %out.MRT_X_72 = ones(dim1,dim2) * NaN;
    %out.MRT_I_72 = ones(dim1,dim2) * NaN;
    out.I_ss_calc = ones(dim1,dim2) * NaN;
    out.X_ss_calc = ones(dim1,dim2) * NaN;
    out.T50_X = ones(dim1,dim2) * NaN;

end

function out = init_out_sens_arrays(out, dim1, dim2)
    % initialize arrays to store 72h values and sensitivities after adjustments
        
    % initialize arrays to store 72h sensitivities after adjustments
    out.sens_X_72 = ones(dim1,dim2) * NaN;
    out.sens_I_72 = ones(dim1,dim2) * NaN;
    out.sens_AUC_X_72 = ones(dim1,dim2) * NaN;
    out.sens_AUC_I_72 = ones(dim1,dim2) * NaN;
    out.sens_AUMC_X_72 = ones(dim1,dim2) * NaN;
    out.sens_AUMC_I_72 = ones(dim1,dim2) * NaN;
    %out.sens_MRT_X_72 = ones(dim1,dim2) * NaN;
    %out.sens_MRT_I_72 = ones(dim1,dim2) * NaN;
    out.sens_I_ss_calc = ones(dim1,dim2) * NaN;
    out.sens_X_ss_calc = ones(dim1,dim2) * NaN;
    out.sens_T50_X = ones(dim1,dim2) * NaN;

end

function out = store_out_val_raw_set(out, i, j, sim_time, sim_y, sp, p, conv_factor_ngml)
    % calculate and store output values and sensitivities

    % store values of I and X at 72h
    out.X_72(i,j) = sim_y(sp.X, end);    % ng/mL
    out.I_72(i,j) = sim_y(sp.I, end);    % #/cell
    
    % store AUC at 72h
    out.AUC_X_72(i,j) = trapz(sim_time, sim_y(sp.X,:)); % (ng * h)/mL
    out.AUC_I_72(i,j) = trapz(sim_time, sim_y(sp.I,:)); % (# * h)/cell

    % store AUMC at 72h
    out.AUMC_X_72(i,j) = trapz(sim_time, sim_time .* sim_y(sp.X,:)); % (ng * h^2)/mL
    out.AUMC_I_72(i,j) = trapz(sim_time, sim_time .* sim_y(sp.I,:)); % (# * h^2)/cell
    
    % store MRT at 72h
    %out.MRT_X_72(i) = out.AUMC_X_72(i) / out.AUC_X_72(i);
    %out.MRT_I_72(i) = out.AUMC_I_72(i) / out.AUC_I_72(i);
    
    % calculate analytic steady state for I, X
    out.I_ss_calc(i,j) = p.alpha / (p.beta + p.gamma);
    out.X_ss_calc(i,j) = p.alpha * p.beta * conv_factor_ngml/ (p.delta * (p.beta + p.gamma));
    
    % T50 for X_ss (time when X is halfway to analytic steady state
    out.T50_X(i,j) = sim_time(dsearchn(sim_y(sp.X,:)', out.X_ss_calc(i)/2));


end

function out = store_out_val_rel_sens_set(out, i, j, frac_change)
    
    % calculate relative sensitivity of X, I, and their AUCs to parameter
    out.sens_X_72(i,j) = ((out.X_72(i,j) - out.X_72_0) / out.X_72_0) / frac_change;
    out.sens_I_72(i,j) = ((out.I_72(i,j) - out.I_72_0) / out.I_72_0) / frac_change;
    out.sens_AUC_X_72(i,j) = ((out.AUC_X_72(i,j) - out.AUC_X_72_0) / out.AUC_X_72_0) / frac_change;
    out.sens_AUC_I_72(i,j) = ((out.AUC_I_72(i,j) - out.AUC_I_72_0) / out.AUC_I_72_0) / frac_change;
    
    % calculate relative sensitivity of AUMC, MRT for X,I
    out.sens_AUMC_X_72(i,j) = ((out.AUMC_X_72(i,j) - out.AUMC_X_72_0) / out.AUMC_X_72_0) / frac_change;
    out.sens_AUMC_I_72(i,j) = ((out.AUMC_I_72(i,j) - out.AUMC_I_72_0) / out.AUMC_I_72_0) / frac_change;
    %out.sens_MRT_X_72(i,j) = ((out.MRT_X_72(i,j) - out.MRT_X_72_0) / out.MRT_X_72_0) / frac_change;
    %out.sens_MRT_I_72(i,j) = ((out.MRT_I_72(i,j) - out.MRT_I_72_0) / out.MRT_I_72_0) / frac_change;
    
    % calculate relative sensitivity of I_ss, X_ss, T50 for X_ss
    out.sens_I_ss_calc(i,j) = ((out.I_ss_calc(i,j) - out.I_ss_calc_0) / out.I_ss_calc_0) / frac_change;
    out.sens_X_ss_calc(i,j) = ((out.X_ss_calc(i,j) - out.X_ss_calc_0) / out.X_ss_calc_0) / frac_change;
    out.sens_T50_X(i,j) = ((out.T50_X(i,j) - out.T50_X_0) / out.T50_X_0) / frac_change;

end

function out = init_out_base_arrays(out, dim1, num_sim_pts)
    % initialize arrays to store baseline values when 
    % testing multiple baselines

    % time courses
    out.X_base_tcs_all = ones(num_sim_pts, dim1, 1) * NaN;
    out.I_base_tcs_all = ones(num_sim_pts, dim1, 1) * NaN;

    % species amounts at 72h
    out.X_72_0 = ones(dim1, 1) * NaN;
    out.I_72_0 = ones(dim1, 1) * NaN; % #/cell
    % species AUC at 72h
    out.AUC_X_72_0 = ones(dim1, 1) * NaN;
    out.AUC_I_72_0 = ones(dim1, 1) * NaN;
    % species AUMC at 72h
    out.AUMC_X_72_0 = ones(dim1, 1) * NaN;
    out.AUMC_I_72_0 = ones(dim1, 1) * NaN;
    % species MRT at 72h - not sure this works without AUC_inf, AUMC_inf
    %out.MRT_X_72_0 = out.AUMC_X_72_0 / out.AUC_X_72_0; % h
    %out.MRT_I_72_0 = out.AUMC_I_72_0 / out.AUC_I_72_0; % h
    % analytic steady state for I, X
    out.I_ss_calc_0 = ones(dim1, 1) * NaN;
    out.X_ss_calc_0 = ones(dim1, 1) * NaN;
    % T50 for X_ss (time when X is halfway to analytic steady state
    out.T50_X_0 = ones(dim1, 1) * NaN;

end

function out = calc_baseline_outs_diffbase(out, h, sim_time, ...
    sim_y, sp, p, conv_factor_ngml)
    % find baseline 72h output values

    % save baseline X,I time courses
    out.X_base_tcs_all(:,h) = sim_y(sp.X,:);
    out.I_base_tcs_all(:,h) = sim_y(sp.I,:);
    
    % species amounts at 72h
    out.X_72_0(h) = sim_y(sp.X, end); % ng/mL
    out.I_72_0(h) = sim_y(sp.I, end); % #/cell
    % species AUC at 72h
    out.AUC_X_72_0(h) = trapz(sim_time, sim_y(sp.X,:)); % (ng * h)/mL
    out.AUC_I_72_0(h) = trapz(sim_time, sim_y(sp.I,:)); % (# * h)/cell
    % species AUMC at 72h
    out.AUMC_X_72_0(h) = trapz(sim_time, sim_time .* sim_y(sp.X,:)); % (ng * h^2)/mL
    out.AUMC_I_72_0(h) = trapz(sim_time, sim_time .* sim_y(sp.I,:)); % (# * h^2)/cell
    % species MRT at 72h - not sure this works without AUC_inf, AUMC_inf
    %out.MRT_X_72_0 = out.AUMC_X_72_0 / out.AUC_X_72_0; % h
    %out.MRT_I_72_0 = out.AUMC_I_72_0 / out.AUC_I_72_0; % h
    % analytic steady state for I, X
    out.I_ss_calc_0(h) = p.alpha / (p.beta + p.gamma);
    out.X_ss_calc_0(h) = p.alpha * p.beta * conv_factor_ngml/ (p.delta * (p.beta + p.gamma));
    % T50 for X_ss (time when X is halfway to analytic steady state
    out.T50_X_0(h) = sim_time(dsearchn(sim_y(sp.X,:)', out.X_ss_calc_0(h)/2));

end

function out = store_out_val_raw_set_diffbase(out, out_index, sim_time, sim_y, sp, p, conv_factor_ngml)
    % calculate and store output values and sensitivities

    % store values of I and X at 72h
    out.X_72(out_index,1) = sim_y(sp.X, end);    % ng/mL
    out.I_72(out_index,1) = sim_y(sp.I, end);    % #/cell
    
    % store AUC at 72h
    out.AUC_X_72(out_index,1) = trapz(sim_time, sim_y(sp.X,:)); % (ng * h)/mL
    out.AUC_I_72(out_index,1) = trapz(sim_time, sim_y(sp.I,:)); % (# * h)/cell

    % store AUMC at 72h
    out.AUMC_X_72(out_index,1) = trapz(sim_time, sim_time .* sim_y(sp.X,:)); % (ng * h^2)/mL
    out.AUMC_I_72(out_index,1) = trapz(sim_time, sim_time .* sim_y(sp.I,:)); % (# * h^2)/cell
    
    % store MRT at 72h
    %out.MRT_X_72(i) = out.AUMC_X_72(i) / out.AUC_X_72(i);
    %out.MRT_I_72(i) = out.AUMC_I_72(i) / out.AUC_I_72(i);
    
    % calculate analytic steady state for I, X
    out.I_ss_calc(out_index,1) = p.alpha / (p.beta + p.gamma);
    out.X_ss_calc(out_index,1) = p.alpha * p.beta * conv_factor_ngml/ (p.delta * (p.beta + p.gamma));
    
    % T50 for X_ss (time when X is halfway to analytic steady state
    out.T50_X(out_index,1) = sim_time(dsearchn(sim_y(sp.X,:)', ...
        out.X_ss_calc(out_index,1)/2));


end

function out = store_out_val_rel_sens_set_diffbase(out, h, out_index, frac_change)


    % calculate relative sensitivity of X, I, and their AUCs to parameter
    out.sens_X_72(out_index,1) = ((out.X_72(out_index,1) - out.X_72_0(h)) / out.X_72_0(h)) / frac_change;
    out.sens_I_72(out_index,1) = ((out.I_72(out_index,1) - out.I_72_0(h)) / out.I_72_0(h)) / frac_change;
    out.sens_AUC_X_72(out_index,1) = ((out.AUC_X_72(out_index,1) - out.AUC_X_72_0(h)) / out.AUC_X_72_0(h)) / frac_change;
    out.sens_AUC_I_72(out_index,1) = ((out.AUC_I_72(out_index,1) - out.AUC_I_72_0(h)) / out.AUC_I_72_0(h)) / frac_change;
    
    % calculate relative sensitivity of AUMC, MRT for X,I
    out.sens_AUMC_X_72(out_index,1) = ((out.AUMC_X_72(out_index,1) - out.AUMC_X_72_0(h)) / out.AUMC_X_72_0(h)) / frac_change;
    out.sens_AUMC_I_72(out_index,1) = ((out.AUMC_I_72(out_index,1) - out.AUMC_I_72_0(h)) / out.AUMC_I_72_0(h)) / frac_change;
    %out.sens_MRT_X_72(out_index,1) = ((out.MRT_X_72(out_index,1) - out.MRT_X_72_0(h)) / out.MRT_X_72_0(h)) / frac_change;
    %out.sens_MRT_I_72(out_index,1) = ((out.MRT_I_72(out_index,1) - out.MRT_I_72_0(h)) / out.MRT_I_72_0(h)) / frac_change;
    
    % calculate relative sensitivity of I_ss, X_ss, T50 for X_ss
    out.sens_I_ss_calc(out_index,1) = ((out.I_ss_calc(out_index,1) - out.I_ss_calc_0(h)) / out.I_ss_calc_0(h)) / frac_change;
    out.sens_X_ss_calc(out_index,1) = ((out.X_ss_calc(out_index,1) - out.X_ss_calc_0(h)) / out.X_ss_calc_0(h)) / frac_change;
    out.sens_T50_X(out_index,1) = ((out.T50_X(out_index,1) - out.T50_X_0(h)) / out.T50_X_0(h)) / frac_change;

end

function [bg_threshold_ub, out, p_mod_all, sim_time, Y_all] = find_bg_threshold(betas, gammas, ...
    bg_threshold_ub, sp, p_mod, prod_case, secr_end_time, time_interval,  ...
    conv_factor_ngml, status_message,  param_names, p_mod_all, Y_all)


        %% Run experimental cases
        for i = 1:length(betas)
            for j = 1:length(gammas)
                
                % update beta, gamma
                p_mod.beta = betas(i);
                p_mod.gamma = gammas(j);

                % calculate index of this parameter set
                out_index = (i-1)*length(gammas) + j;

                % store p_mod values for data analysis
                for x = 1:length(param_names)
                    p_mod_all(x, out_index) = p_mod.(param_names{x});
                end

                % Simulation to steady state
                try
                    [sim_time, sim_y] = sim_secr_dde_to_ss(sp, p_mod, prod_case, secr_end_time,...
                        time_interval, conv_factor_ngml, status_message);

                    % store phenomenological parameters (72h values, AUCs, etc)
                    % out = store_out_val_raw_set(out, i, j, sim_time, sim_y, sp, p, conv_factor_ngml, factors_to_try(j));
                
                    % store values of I and X at 72h
                    out.X_ss(i,j) = sim_y(sp.X, end);    % ng/mL
                    out.I_ss(i,j) = sim_y(sp.I, end);    % #/cell

                    % if the simulation failed because I is negative
                    if sim_y(sp.I, end) < 0

                        % calculate beta + gamma
                        bg = betas(i) + gammas(j);

                        % initialize b+g threshold if needed
                        if isnan(bg_threshold_ub)
                            bg_threshold_ub = bg;
                        end
                        
                        % update b+g threshold if needed
                        if bg < bg_threshold_ub
                            bg_threshold_ub = bg;
                        end
                    end

                catch

                    fprintf("Failed parameter sets: %d x %d \n", betas(i), gammas(j))
                    continue
                
                end

            end
        end

        %sim_time_out = sim_time;

end