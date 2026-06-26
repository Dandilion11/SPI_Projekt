%Runs the optimization algorithm (e.g., lsqnonlin or fmincon) minimizing the weighted least squares objective function (Task 4).
% main_04_parameter_fitting.m
clear; clc;

%% 1. Load Preprocessed Data
load('Processed_Batch_Data.mat');
t_messung = TrainData.Biomasse.t;
y_bio = TrainData.Biomasse.y;
y_glc = TrainData.Glucose.y;
var_bio = TrainData.Biomasse.var;
var_glc = TrainData.Glucose.var;

%% 2. System Constants & Initial Conditions
V = 10; % Constant volume for simple batch (verify exact value with reactor data)

% Load dynamic volume
%V = TrainData.V; 


m_X0 = y_bio(1) * V;
m_Glc0 = y_glc(1) * V;
x0 = [m_X0; m_Glc0];

%% 3. Optimization Setup
% Parameter vector: [mu_max, K_Glc, Y_Glc_X]
p0 = [0.3, 0.5, 0.15]; 
lb = [0.01, 0.01, 0.01]; 
ub = [1.0, 5.0, 1.0];

% Select Kinetic Model (1 = Monod, 2 = Tessier, 3 = Moser)
model_type = 1; 

%% 4. Define Objective Function (Anonymous Function Wrapper)
% Wraps the custom WLS error calculation to only expose 'p' to fmincon
obj_fun = @(p) calculate_wls_error(p, x0, t_messung, y_bio, y_glc, var_bio, var_glc, model_type, V);

%% 5. Execute Optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[p_opt, fval] = fmincon(obj_fun, p0, [], [], [], [], lb, ub, [], options);

%% 6. Output Results
fprintf('\n--- Optimization Complete ---\n');
fprintf('Final Weighted Least Squares Error: %.4f\n', fval);
fprintf('mu_max  = %.4f h^-1\n', p_opt(1));
fprintf('K_Glc   = %.4f g/L\n', p_opt(2));
fprintf('Y_Glc_X = %.4f g/g\n', p_opt(3));

%% 7. Visual Validation (Dark Theme Fix)
% Generate a continuous time vector for a smooth simulation line
t_sim = linspace(0, 20, 200);

% Simulate the system using the optimized parameters (p_opt)
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
[~, X_opt] = ode15s(@(t, x) ode_task2a_batch(t, x, p_opt, model_type, V), t_sim, x0, options);

% Convert simulated mass back to concentration
c_X_sim   = X_opt(:, 1) / V;
c_Glc_sim = X_opt(:, 2) / V;

% Plot Results
figure('Name', 'Task 2a: Parameter Identification (Monod)', 'Position', [150, 150, 900, 600]);

% Biomass Plot
subplot(2, 1, 1);
% 'wo' (weiß) anstelle von 'ko' (schwarz), und dicke Cyan-Linie für Simulation
errorbar(t_messung, y_bio, sqrt(var_bio), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_X_sim, 'c-', 'LineWidth', 2.5);
title(sprintf('Biomass Fit (Error: %.1f)', fval), 'Color', 'w');
xlabel('BatchAge (h)', 'Color', 'w');
ylabel('c_X (g/L)', 'Color', 'w');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w'); % Dunkler Hintergrund
legend('Experimental Data \pm \sigma', 'Optimized Model', 'Location', 'NorthWest', 'TextColor', 'w');
grid on;

% Glucose Plot
subplot(2, 1, 2);
% 'wo' (weiß), und dicke gelbe Linie für Simulation
errorbar(t_messung, y_glc, sqrt(var_glc), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_Glc_sim, 'y-', 'LineWidth', 2.5);
title('Glucose Fit', 'Color', 'w');
xlabel('BatchAge (h)', 'Color', 'w');
ylabel('c_{Glc} (g/L)', 'Color', 'w');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w'); % Dunkler Hintergrund
legend('Experimental Data \pm \sigma', 'Optimized Model', 'Location', 'NorthEast', 'TextColor', 'w');
grid on;

sgtitle(sprintf('Task 2a: Model Type %d Optimization Results', model_type), 'Color', 'w');