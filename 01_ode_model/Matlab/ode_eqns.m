function dydt = ode_eqns(t, y, sp, p, prod_case, conv_factor)

% SIMPLE_EQNS   Equations file for ode15s simulating sFlt1 trafficking with
%               4 equations in one cell type
%
% EXAMPLE USAGE (within ode15s)
% -----------------------------
% [T_pulse, Y_pulse] = ode15s(@simple_eqns, pulse_tspan, y0, ...
%                             ode_options, sp, p, prod_case, conv_factor);
%
% OUTPUT (to ode15s)
% ------------------
% dydt:     rates of change of concentrations of molecules at time t
%
% INPUT (from ode15s)
% -------------------
% t:            current time
% y:            current molecule amounts/concentrations (see units below)
% sp:           structure mapping molecule types to index numbers
% p:            structure containing rate constant parameter values
% prod_case:    string listing production case for assigning transfer function
% conv_factor:  conversion factor for intracellular -> extracellular rxns
%
% UNITS
% -----
% Intracellular:    #/cell
% Extracellular:    #/cell
%

% Initialize dydt
dydt = zeros(length(fieldnames(sp)), 1);

% Assign transfer function value
switch prod_case
    case "on"    % max production
        trans_fun = 1;
    case "off"   % zero production
        trans_fun = 0;
end

%% EQUATIONS

% intracellular sFlt1
dydt(sp.I) = p.alpha * trans_fun - p.beta * y(sp.I) - p.gamma * y(sp.I);

% intracellularly degraded sFlt1
dydt(sp.I_deg) = p.gamma * y(sp.I);

% extracellular sFlt1
dydt(sp.X) = p.beta * y(sp.I) * conv_factor - p.delta * y(sp.X);

% extracellularly degraded sFlt1
dydt(sp.X_deg) = p.delta * y(sp.X);