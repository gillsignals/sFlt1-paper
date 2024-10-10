function [T_final, Y_final, X_frac, I_fc] = sim_secr_ode(sp, p, ...
    secr_end_time, time_interval, conv_factor)

% SIM_SECR_ODE  Simulate response to media change for sFlt1 trafficking 
%               with specified number of equations in one cell type.
%   
% INPUT
% -----
%
% prod_decay:   specifies the rate of decay of labeled protein production 
%               as a function of time since the start of the chase.
%               Available options are:
%
%       "on":          production fully off during chase (Heaviside function)
%       "chase":        production exponentially decays during chase
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

% parameters for finding steady state
t_int = 20;                 % interval of time to cover per iteration
count = 1;                  % used to track iterations before steady state found
tspan = 0:time_interval:t_int;          % specify timespan for simulation (minutes)
max_delta = Inf;            % initialize max_delta to Inf (will reduce as simulation converges to steady state)
                
% set ODE solver options
ode_options = odeset('AbsTol', 1e-12, ...
    'RelTol', 1e-8, ...
    'NonNegative', 1:length(fieldnames(sp)), ...
    'InitialStep', 1e-2);

% set all initial species values to 0
y0 = zeros(length(fieldnames(sp)), 1);

% specify production fully on for entire simulation
prod_case = "on";

% Run ode15s to generate simulation results using associated equations file
[T, Y] = ode15s(@ode_eqns, tspan, y0, ode_options, sp, p, prod_case, conv_factor);
                
while max_delta > .005

    % update simulation conditions
    count = count + 1;                              % increment count
    tspan = (t_int * (count - 1)):time_interval:(t_int * count);   % calculate new time range
    y0_new = Y(end,:);

    % if count is too high, say there is no steady state
    if count > 100

        disp("Infinite sadness =(: Steady state can't be reached.")

        infinite_sadness = true;

       return

    end

    % simulate an additional time interval

    % Run ode15s to generate simulation results using associated equations file
    [~, Y_new] = ode15s(@ode_eqns, tspan, y0_new, ode_options, ...
        sp, p, prod_case, conv_factor);

    % add simulated data from new interval to final output
    T = [T; tspan(2:end)'];
    Y = [Y; Y_new];

    % specify species to check for steady state (not degraded)
    % UPDATE 1/26/23 - only check Si, not Sx
    ind = [sp.I];

    % calculate largest error in concentration for non-degraded R1 species
    current_concs = Y_new(end, ind);   
    prior_concs = y0_new(ind);
    delta = abs((current_concs - prior_concs) ./ prior_concs);
    max_delta = max(delta(delta < Inf));

end

% set new history / initial conditions to steady state
y0_ss = Y(end,:);

% find baseline steady state values
X_ss = Y(end, sp.X);
I_ss = Y(end, sp.I);

% simulate media change by setting X to 0
y0_ss(sp.X) = 0;

% reset degraded species too
y0_ss(sp.X_deg) = 0;
y0_ss(sp.I_deg) = 0;

% update simulation time to experimental time
T_final = 0:time_interval:secr_end_time;

% simulate time interval 
    [~, Y_final] = ode15s(@ode_eqns, T_final, y0_ss, ode_options, ...
        sp, p, prod_case, conv_factor);

% convert X to fraction of 24h signal
X_24h = Y_final(dsearchn(T_final', 24), sp.X);
X_frac = Y_final(:,sp.X) ./ X_24h;

% convert I to fold change from time 0
I_fc = Y_final(:,sp.I) ./ I_ss;


end