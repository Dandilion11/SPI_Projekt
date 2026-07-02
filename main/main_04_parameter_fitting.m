
% main_04_parameter_fitting.m
clear; clc; close all;
%{
%% 1. Load Preprocessed Data
load('Processed_Batch_Data.mat');
t_messung = TrainData.Biomasse.t;
y_bio = TrainData.Biomasse.y;
y_glc = TrainData.Glucose.y;
var_bio = TrainData.Biomasse.var;
var_glc = TrainData.Glucose.var;
%}

% Exp3 = TrainData
% Exp4 = ValData
% mit ValData läuft die PI
% Ich glaube bei Exp3 wird Glucose zugeführt (ca. stunde 8). Dadurch wäre
% das kein geeigenetes experiment mehr für das simpe model.
% Wenn tatsächlich glucose gefüttert wird, könnten wir bei exp3 die daten
% bis zum füttern nehmen


%% 1. Load Preprocessed Data
projectRoot = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(projectRoot,'Daten/Daten_Processed/Processed_Batch_Data.mat'));
addpath(fullfile(projectRoot,'Modelle'),'-begin');

% SWAPPED: Using ValData (Experiment 04) for training
t_messung = ValData.Biomasse.t;
y_bio = ValData.Biomasse.y;
y_glc = ValData.Glucose.y;
var_bio = ValData.Biomasse.var;
var_glc = ValData.Glucose.var;

%% 2. System Constants & Initial Conditions
% Modell1 expects concentrations (g/L), volume V is removed.
x0 = [y_bio(1); y_glc(1)];

%% 3. Optimization Setup
% Parameter vector: [mumax, KS, YXS]
p0 = [0.3, 0.5, 0.15]; 
lb = [0.01, 0.01, 0.01]; 
ub = [1.0, 5.0, 1.0]; %KS upper bound

% Flags for Modell1
kinetic = 3; % 3 = Monod
withOxygen = false;

%% 4. Define Objective Function
obj_fun = @(p) calculate_wls_error(p, x0, t_messung, y_bio, y_glc, var_bio, var_glc, kinetic, withOxygen);

%% 5. Execute Optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[p_opt, fval] = fmincon(obj_fun, p0, [], [], [], [], lb, ub, [], options);

%% 6. Output Results
fprintf('\n--- Optimization Complete ---\n');
fprintf('Final WLS Error: %.4f\n', fval);
fprintf('mumax  = %.4f h^-1\n', p_opt(1));
fprintf('KS     = %.4f g/L\n', p_opt(2));
fprintf('YXS    = %.4f g/g\n', p_opt(3));

scriptDir = fileparts(mfilename('fullpath'));
saveDir = fullfile(scriptDir, '..', 'Daten');
savePath = fullfile(saveDir, 'p_opt.mat');
save(savePath,"p_opt");

%% 7. Visual Validation
t_sim = linspace(0, 20, 200);
options_ode = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
[~, X_opt] = ode15s(@(t, x) Modell1(t, x, p_opt, kinetic, withOxygen), t_sim, x0, options_ode);

c_X_sim   = X_opt(:, 1);
c_Glc_sim = X_opt(:, 2);

figure('Name', 'Task 2a: Parameter Identification (Modell1)', 'Position', [150, 150, 900, 600]);

subplot(2, 1, 1);
errorbar(t_messung, y_bio, sqrt(var_bio), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_X_sim, 'c-', 'LineWidth', 2.5);
title(sprintf('Biomass Fit (Error: %.1f)', fval), 'Color', 'w');
ylabel('c_X (g/L)', 'Color', 'w');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w');
legend('Experimental Data \pm \sigma', 'Optimized Model', 'TextColor', 'w');
grid on;

subplot(2, 1, 2);
errorbar(t_messung, y_glc, sqrt(var_glc), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_Glc_sim, 'y-', 'LineWidth', 2.5);
title('Glucose Fit', 'Color', 'w');
xlabel('BatchAge (h)', 'Color', 'w');
ylabel('c_{Glc} (g/L)', 'Color', 'w');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w');
legend('Experimental Data \pm \sigma', 'Optimized Model', 'TextColor', 'w');
grid on;

%% Local Objective Function
function J = calculate_wls_error(p, x0, t_mes, y_bio, y_glc, var_bio, var_glc, kinetic, withOxygen)
    % Force ode15s to evaluate exactly at measurement timestamps
    try
        [~, X_sim] = ode15s(@(t, x) Modell1(t, x, p, kinetic, withOxygen), t_mes, x0);
        c_X_sim = X_sim(:,1);
        c_Glc_sim = X_sim(:,2);
        
        err_bio = sum(((y_bio - c_X_sim).^2) ./ var_bio);
        err_glc = sum(((y_glc - c_Glc_sim).^2) ./ var_glc);
        
        J = err_bio + err_glc;
    catch
        J = 1e6; % Penalty for integration failure
    end
end