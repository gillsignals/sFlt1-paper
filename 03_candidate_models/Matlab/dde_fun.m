function dydt = dde_fun(t, y, Z, sp, p, prod_case, mat_case, internalize_case, conv_factor)
    
% OUTPUT
% dydt: rates of change of concentrations of molecules at time t

% INPUT
% t: current time
% y: current molecule amounts/concentrations (see units below)
% Z: function history - inside dde23, z(:,j) is molecule concs at time t-z(j)
%    example: if lags = (2, 4), then z(:,1) means "conc at time t-2" and
%       z(:,2) means "conc at time t-4"
%    works if you have multiple species: z(1,1) is "conc of y_1 at t-lag_1",
%    z(i,j) is "conc of y_i at t-lag_j 
% sp: structure mapping molecule types to index numbers
% p: structure containing rate constant parameter values
% prod_case: string listing production case for assigning transfer function
%   "off":          production fully off during chase (Heaviside function)
%   "lin_decay":    production linearly decays over time during chase
%   "exp_decay":    production exponentially decays during chase
%   "delay_off":    production fully off after a cutoff time during chase
%                   (delayed Heaviside function)
% num_cells: number of cells (to convert intra/extracellular amounts
% conv_factor: conversion factor for intracellular -> extracellular rxns

% UNITS:
% Intracellular: #/cell
% Extracellular: ng/mL

% find y(t - ylag)
ylag = Z(:,1);

% Assign transfer function value for production reaction
switch prod_case
    case "on"    % max production
        trans_fun = 1;
    case "off"   % zero production
        trans_fun = 0;
    case "lin_decay" % linearly decaying production
        trans_fun = max(1 - p.kappa * t, 0); % 1 - kt    
    case "exp_decay" % exponentially decaying production
        trans_fun = exp(1) ^ -(p.kappa * t); % e^-kt
    case "delay_off" % zero production after time = p.kappa
        if t < p.kappa
            trans_fun = 1;
        else
            trans_fun = 0;
        end
end

% specify production flux based on prod_case (#/cell/h)
prod_flux = p.alpha * trans_fun;

% specify secretion flux and intracellular degradation flux based on mat_case (#/cell/h)
switch mat_case
    case 'no'
        secr_flux = p.beta * y(sp.I);
        intdeg_flux = p.gamma * y(sp.I);
    case 'yes'
        secr_flux = p.beta * ylag(sp.I);
        intdeg_flux = p.gamma * ylag(sp.I);
end

% specify internalization flux based on internalize_case (#/cell/h)
switch internalize_case
    case 'no'
        internalize_flux = 0;
    case 'yes'
        internalize_flux = p.epsilon * y(sp.X) / conv_factor;
end

% specify external degradation flux (#/cell/h)
extdeg_flux = p.delta * y(sp.X) / conv_factor;

% Initialize dydt
dydt = zeros(1, length(fieldnames(sp)));
    
% I change = production - delayed secretion - delayed degradation +
% internalization (#/cell/h)
if y(sp.I) >= 0
    dydt(sp.I) = prod_flux - secr_flux - intdeg_flux + internalize_flux;
else
    dydt(sp.I) = 0;
end

% degraded I change = delayed degradation 
dydt(sp.I_deg) = intdeg_flux;
    
% X change = delayed secretion - degradation - internalization (ng/mL/h)
dydt(sp.X) = (secr_flux - extdeg_flux - internalize_flux) * conv_factor;
    
% degraded X change = degradation (#/cell/h)
dydt(sp.X_deg) = extdeg_flux;

% convert to column vector (gets auto-converted to row vector after assignment)
dydt = dydt';

end