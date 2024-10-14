function [T_full, Y_full, X_frac, I_fc] = sim_pulse_chase_dde_v2(prod_decay, mat_case, ...
        internalize_case, sp, p, pulse_length, chase_length, time_interval, conv_factor);

% SIM_PULSE_CHASE   Simulate a pulse-chase experiment for sFlt1 trafficking 
%                   with specified number of equations in one cell type.
%   
% INPUT
% -----
%
% prod_decay:   specifies the rate of decay of labeled protein production 
%               as a function of time since the start of the chase.
%               Available options are:
%
%       "off":          production fully off during chase (Heaviside function)
%       "lin_decay":    production linearly decays over time during chase
%       "exp_decay":    production exponentially decays during chase
%       "delay_off":    production fully off after a cutoff time during chase
%                       (delayed Heaviside function)
%
% pulse_length:     length of pulse (in hours)
%
% chase_length:     length of chase (in hours)
%
% time_interval:    time point intervals at which to return species values;
%                   required to make output vectors a predetermined length
% 
% conv_factor:      conversion factor for #/cell to pmol
%
% See also DDE_FUN, DDESET, DDE23.


%% 

% set ODE solver options
dde_options = ddeset('AbsTol', 1e-6, ...
    'RelTol', 1e-3, ...
    'InitialStep', 1e-2);
 %   'NonNegative', 1:length(fieldnames(sp)), ... % not available for dde?


% set all initial species values to 0
y0 = zeros(length(fieldnames(sp)), 1);

% specify time spans for pulse and chase
pulse_tspan = -pulse_length:time_interval:0; % 20 minutes
chase_tspan = 0:time_interval:chase_length;    % 10 hours

%% SIMULATE PULSE

% specify production fully on for pulse
prod_case = "on";

history_0 = zeros(1, length(fieldnames(sp)));

% starting from all species at 0, simulate 20 minutes of production

sol_pulse = dde23(@dde_fun_v2, p.tau, history_0, pulse_tspan, dde_options, ...
    sp, p, prod_case, mat_case, internalize_case, conv_factor);

% extract time and species values from pulse solution
Y_pulse = deval(sol_pulse, pulse_tspan);

%fig = figure("Name","Temp pulse plot");

% Extracellular sFlt1 (Sx) amount
%subplot(2,1,1);
%plot(T_pulse, Y_pulse(sp.X,:), 'LineWidth', 3)
%title("Extracellular S Amount (S_x)")
%xlabel("Time (h)")
%ylabel("Amount (pmol)")

% Intracellular sFlt1 (Si) amount
%subplot(2,1,2)
%plot(T_pulse, Y_pulse(sp.I, :), 'LineWidth', 3)
%title("Intracellular S Amount (S_i)")
%xlabel("Time (h)")
%ylabel("Amount (#/cell)")

%% SIMULATE CHASE

% update production case depending on specified decay type
prod_case = prod_decay;

% set initial conditions for chase to conditions at end of pulse
y0_chase = sol_pulse.y(:,end);
        
% simulate media change by setting Sx to 0
y0_chase(sp.X) = 0;

% update DDE options to note locations of discontinuities
dde_options = ddeset(dde_options, 'InitialY', y0_chase);

% starting after pulse, simulate 10 hours with no production
sol_chase = dde23(@dde_fun_v2, p.tau, sol_pulse, chase_tspan, dde_options, ...
    sp, p, prod_case, mat_case, internalize_case, conv_factor);

% extract time and species values from chase solution
% note the eps() - limit from the right for the discontinuity at 0
Y_chase = deval(sol_chase, [eps, chase_tspan(2:end)]);

%% MERGE PULSE-CHASE TIME COURSES

% combine times, treating start of chase as t=0
% include two time points before pulse to distinguish pulse effects
T_full = [-pulse_length*3, -pulse_length, pulse_tspan, chase_tspan];
        
% combine Y values
% including two zeros at times before pulse to distinguish pulse effects
Y_full = [zeros(length(fieldnames(sp)),2), Y_pulse, Y_chase];

% convert X to fraction of 8h signal
X_8h = Y_full(sp.X, dsearchn(T_full', 8));
X_frac = Y_full(sp.X,:) ./ X_8h;

% convert I to fold change from end of pulse
I_fc = Y_full(sp.I,:) ./ y0_chase(sp.I);


end