function dxdt = ode_task2a_batch(t, x, params, model_type, V)
% ODE system for Task 2a: Simple Batch Phase
%
% Inputs:
% t          : Time (h) - required by ODE solver
% x          : State vector [m_X; m_Glc] in (g)
% params     : Parameter vector 
%              [mu_max, K_Glc, Y_Glc_X, n]
% model_type : Integer (1 = Monod, 2 = Tessier, 3 = Moser)
% V          : Constant reactor volume (L)
%
% Outputs:
% dxdt       : Derivatives [dm_X/dt; dm_Glc/dt]

% 1. Unpack States
m_X   = x(1);
m_Glc = x(2);

% 2. Calculate Concentrations (c = m / V)
c_Glc = m_Glc / V;

% 3. Unpack Parameters
mu_max  = params(1);
K_Glc   = params(2);
Y_Glc_X = params(3);

% 3. Unpack Parameters
mu_max  = params(1);
K_Glc   = params(2);
Y_Glc_X = params(3);

% Safely extract 'n' only if provided (for Moser)
if length(params) == 4
    n = params(4);
else
    n = 1; % Fallback for Monod and Tessier
end

% 4. Calculate Specific Growth Rate (r_X) based on selected model
switch model_type
    case 1
        % Monod Kinetics
        r_X = mu_max * (c_Glc / (c_Glc + K_Glc));

    case 2
        % Tessier Kinetics
        r_X = mu_max * (1 - exp(-c_Glc / K_Glc));

    case 3
        % Moser Kinetics
        r_X = mu_max * (c_Glc^n / (c_Glc^n + K_Glc));

    otherwise
        error('Invalid model_type. Use 1 (Monod), 2 (Tessier), or 3 (Moser).');
end

% 5. Formulate Absolute Mass Balances
dm_X_dt   = r_X * m_X;
dm_Glc_dt = -(1 / Y_Glc_X) * r_X * m_X;

% 6. Return Derivative Vector
dxdt = [dm_X_dt; dm_Glc_dt];
end