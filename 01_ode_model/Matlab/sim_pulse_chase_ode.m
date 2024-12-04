function [T_full, Y_full, X_frac, I_fc] = sim_pulse_chase_ode(sp, p, ...
    pulse_length, chase_length, time_interval, conv_factor)

% SIM_PULSE_CHASE_ODE   Simulate a pulse-chase experiment for sFlt1 trafficking 
%                       with specified number of equations in one cell type.
%   
% INPUT
% -----
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
% See also ODE_EQNS, ODE15S.


%% 

% set ODE solver options
ode_options = odeset('AbsTol', 1e-12, ...
    'RelTol', 1e-8, ...
    'NonNegative', 1:length(fieldnames(sp)), ...
    'InitialStep', 1e-2);

% set all initial species values to 0
y0 = zeros(length(fieldnames(sp)), 1);

% specify time spans for pulse and chase
pulse_tspan = -pulse_length:time_interval:0; % 20 minutes
chase_tspan = 0:time_interval:chase_length;    % 10 hours

%% SIMULATE PULSE

% specify production fully on for pulse
prod_case = "on";

% starting from all species at 0, simulate 20 minutes of production
% Run ode15s to generate simulation results using associated equations file
[~, Y_pulse] = ode15s(@ode_eqns, pulse_tspan, y0, ...
    ode_options, sp, p, prod_case, conv_factor);

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

% turn off production
prod_case = "off";

% set initial conditions for chase to conditions at end of pulse
y0_chase = Y_pulse(end, :);
        
% simulate media change by setting X to 0
y0_chase(sp.X) = 0;

[~, Y_chase] = ode15s(@ode_eqns, chase_tspan, y0_chase, ...
            ode_options, sp, p, prod_case, conv_factor);

%% MERGE PULSE-CHASE TIME COURSES

% combine times, treating start of chase as t=0
% include two time points before pulse to distinguish pulse effects
T_full = [-pulse_length*3, -pulse_length, pulse_tspan, chase_tspan];
        
% combine Y values
% including two zeros at times before pulse to distinguish pulse effects
Y_full = [zeros(2,length(fieldnames(sp))); Y_pulse; Y_chase];

% convert X to fraction of 8h signal
X_8h = Y_full(dsearchn(T_full', 8), sp.X);
X_frac = Y_full(:,sp.X) ./ X_8h;

% convert I to fold change from end of pulse
I_fc = Y_full(:,sp.I) ./ y0_chase(sp.I);

end