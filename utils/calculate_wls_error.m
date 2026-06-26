function J = calculate_wls_error(p, x0, t_messung, y_bio_mess, y_glc_mess, var_bio, var_glc, model_type, V)
% Objective function for Weighted Least Squares (WLS) parameter fitting
%
% Inputs:
% p          : Current parameter guess [mu_max, K_Glc, Y_Glc_X]
% x0         : Initial states (mass) [m_X0; m_Glc0]
% t_messung  : Time vector of the offline measurements
% y_*_mess   : Measured concentration data (g/L)
% var_* : Calculated variance for weighting
% model_type : 1=Monod, 2=Tessier, 3=Moser
% V          : Reactor volume (L)

% 1. Setup ODE Solver Options
% Strict tolerances prevent the solver from making integration errors
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

% 2. Run the Simulation
% By passing t_messung directly as the time span, ode15s forces the 
% simulation to output data at the exact timestamps of your measurements.
try
    [~, X_sim] = ode15s(@(t, x) ode_task2a_batch(t, x, p, model_type, V), t_messung, x0, options);
catch
    % If the optimizer guesses wild parameters that crash the ODE solver,
    % return a massive error to tell fmincon "bad guess, go the other way."
    J = 1e6; 
    return;
end

% 3. Extract simulated mass and convert back to concentration (c = m/V)
% X_sim(:, 1) is all rows of the first column (Biomass mass)
c_X_sim   = X_sim(:, 1) / V; 
c_Glc_sim = X_sim(:, 2) / V;

% Ensure all data vectors are column vectors (prevents matrix dimension errors)
y_bio_mess = y_bio_mess(:);
y_glc_mess = y_glc_mess(:);
var_bio    = var_bio(:);
var_glc    = var_glc(:);
c_X_sim    = c_X_sim(:);
c_Glc_sim  = c_Glc_sim(:);

% 4. Calculate the Error (Weighted Least Squares)
% Formula: Sum of ((Simulated - Measured)^2 / Variance)
error_bio = sum(((c_X_sim - y_bio_mess).^2) ./ var_bio);
error_glc = sum(((c_Glc_sim - y_glc_mess).^2) ./ var_glc);

% 5. Total Objective Value (J)
J = error_bio + error_glc;

% Safety catch for mathematically invalid numbers (NaN or Infinity)
if isnan(J) || isinf(J)
    J = 1e6;
end
end