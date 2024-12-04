function [T_norm, Y_full, X_norm, I_norm] = sim_on_off_dde(sp, p, ...
    pulse_length, chase_length, time_interval, t50_X, I_ss, X_ss, conv_factor)

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

% set DDE solver options
dde_options = ddeset('AbsTol', 1e-4, ...
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

% define DDE history
history_0 = zeros(1, length(fieldnames(sp)));

% starting from all species at 0, simulate a period of production
% Run ode15s to generate simulation results using associated equations file

sol_on = dde23(@dde_eqns, p.tau, history_0, pulse_tspan, dde_options, ...
    sp, p, prod_case, conv_factor);

Y_pulse = deval(sol_on, pulse_tspan);

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
y0_chase = Y_pulse(:, end);

% run simulation
sol_off = dde23(@dde_eqns, p.tau, y0_chase', chase_tspan, dde_options, ...
    sp, p, prod_case, conv_factor);

% extract time and species values from simulation after media change
Y_chase = deval(sol_off, chase_tspan);

% combine times, treating start of chase as t=0
% include two time points before pulse to distinguish pulse effects
T_norm = [-pulse_length*3, -pulse_length, pulse_tspan, chase_tspan] ./ t50_X;
        
% combine Y values
% including two zeros at times before pulse to distinguish pulse effects
Y_full = [zeros(length(fieldnames(sp)),2), Y_pulse, Y_chase];

%% NORMALIZE X,I

I_ss = p.alpha/(p.beta + p.gamma);
X_ss = p.beta * I_ss / p.delta;

% convert X to fraction of steady state 
X_norm = Y_full(sp.X,:) ./ (X_ss .* conv_factor);

% convert I to fraction of steady state
I_norm = Y_full(sp.I,:) ./ I_ss;

end