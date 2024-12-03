function [T_final, Y_final] = sim_secr_dde(sp, p, prod_case, mat_case, internalize_case, ...
    secr_end_time, time_interval, conv_factor, status_message)

% SIM_SECR_DDE  Simulate response to media change for sFlt1 trafficking 
%               with specified number of equations in one cell type.
%   
% INPUT
% -----
%
% prod_case:   specifies the rate of decay of labeled protein production 
%               as a function of time since the start of the chase.
%               Available options are:
%
%       "on":           constitutive production
%       "off":          production fully off during chase (Heaviside function)
%       "lin_decay":    production linearly decays over time during chase
%       "exp_decay":    production exponentially decays during chase
%       "delay_off":    production fully off after a cutoff time during chase
%                       (delayed Heaviside function)
%
% secr_end_time:    time of follow-up for secretion after media change (hours)
%
% time_interval:    time point intervals at which to return species values;
%                   required to make output vectors a predetermined length
% 
% conv_factor:      conversion factor for #/cell to ng/mL
%
% See also DDE_FUN, DDESET, DDE23.


%% FIND THE STEADY STATE OF THE SYSTEM GIVEN CURRENT PARAMETERS

if ismember(status_message, [{'all'}])
    disp("Finding steady state...");
end

% parameters for finding steady state
t_int = 20;                 % interval of time to cover per iteration
count = 1;                  % used to track iterations before steady state found
tspan = 0:time_interval:t_int;          % specify timespan for simulation (minutes)
max_delta = Inf;            % initialize max_delta to Inf (will reduce as simulation converges to steady state)
                
% set DDE solver options
dde_options = ddeset('AbsTol', 1e-4, ...
    'RelTol', 1e-3, ...
    'InitialStep', 1e-2);
 %   'NonNegative', 1:length(fieldnames(sp)), ... % not available for dde?


% set all initial species values to 0
y0 = zeros(length(fieldnames(sp)), 1);

% specify production fully on for entire simulation
%prod_case = "on";

% define DDE history
history_0 = zeros(1, length(fieldnames(sp)));

% Run ode15s to generate simulation results using main_ode function and associated equations file
if ismember(status_message, [{'all'}])
    fprintf('Running simulation from hours %d to %d...\n', tspan(1), tspan(end));
end

sol_find_ss = dde23(@dde_eqns, p.tau, history_0, tspan, dde_options, ...
    sp, p, prod_case, conv_factor);

T = tspan;
Y = deval(sol_find_ss, tspan);
                
while max_delta > .005
    if ismember(status_message, [{'all'}])
        fprintf('Maximum error: %d \n', max_delta);
    end

    % update simulation conditions
    count = count + 1;                              % increment count
    tspan = (t_int * (count - 1)):time_interval:(t_int * count);   % calculate new time range
    y0_new = sol_find_ss.y(:,end);

    % if count is too high, say there is no steady state
    if count > 100

        disp("Infinite sadness =(: Steady state can't be reached.")

        infinite_sadness = true;

       return

    end

    % simulate an additional time interval
    if ismember(status_message, [{'all'}])
        fprintf('Running simulation from hours %d to %d...\n', tspan(1), tspan(end));
    end

    sol_find_ss = dde23(@dde_eqns, p.tau, sol_find_ss, tspan, dde_options, ...
    sp, p, prod_case, conv_factor);

    Y_new = deval(sol_find_ss, tspan(2:end));

    % add simulated data from new interval to final output
    T = [T, tspan(2:end)];
    Y = [Y, Y_new];

    % specify species to check for steady state (not degraded)
    % UPDATE 1/26/23 - only check I, not X (since we're doing media change anyway)
    ind = [sp.I];

    % calculate largest change in concentration for non-degraded species
    current_concs = Y_new(ind, end);   
    prior_concs = y0_new(ind);
    delta = abs((current_concs - prior_concs) ./ prior_concs);
    max_delta = max(delta(delta < Inf));

end
                
if ismember(status_message, [{'all'}])
    disp("Steady state reached!")
end

% set new history / initial conditions to steady state
y0_ss = Y(:, end);

% find baseline steady state values
X_ss = Y(sp.X, end);
I_ss = Y(sp.I, end);
    
if ismember(status_message, [{'all'}])
    disp("Steady state concentrations before perturbation:")
    disp("X_ss = " + y0_ss(sp.X) + " ng/ml)")
    disp("I_ss = " + y0_ss(sp.I) + " #/cell")
end       

% simulate media change by setting X to 0
y0_ss(sp.X) = 0;

% reset degraded species too
y0_ss(sp.X_deg) = 0;
y0_ss(sp.I_deg) = 0;

% update simulation time to experimental time
T_final = 0:time_interval:secr_end_time;

% simulate time interval 
sol = dde23(@dde_eqns, p.tau, y0_ss', T_final, dde_options,...
    sp, p, prod_case, conv_factor);

% extract time and species values from simulation after media change
Y_final = deval(sol, T_final);

end