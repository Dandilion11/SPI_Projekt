% main_02_simulation_test.m
clear; clc; close all;

%% 1. Load Preprocessed Data
scriptDir = pwd;
loadPath = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');
projectRoot = pwd;
addpath(fullfile(projectRoot, '..','Modelle'),'-begin');

rehash;
clear Modell1

if ~isfile(loadPath)
    error('Data not found.');
end
load(loadPath);

t_messung = TrainData.Biomasse.t;
y_bio     = TrainData.Biomasse.y;
y_glc     = TrainData.Glucose.y;

%% 2. Define Initial Conditions and Parameters
% Modell1 expects concentrations (cX, cGlc), not masses.
x0 = [y_bio(1); y_glc(1)];

% Parameters for Monod: [mumax, KS, YXS]
p_guess = [0.3, 0.5, 0.15]; 

% Model Configuration Flags
kinetic = 3; % 3 = Monod
withOxygen = false; % Exclude O2 state and parameters

%% 3. Execute Simulation
t_sim = linspace(0, 20, 200);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

% Call Modell1 with the required flags
[~, X_sim] = ode15s(@(t, x) Modell1(t, x, p_guess, kinetic, withOxygen), t_sim, x0, options);

%% 4. Extract Data
% Division by V is no longer necessary as Modell1 outputs concentrations directly.
c_X_sim   = X_sim(:, 1);
c_Glc_sim = X_sim(:, 2);

%% 5. Plot Results
figure('Name', 'Simulation Test (Modell1)', 'Position', [200, 200, 900, 600]);

subplot(2, 1, 1);
plot(t_messung, y_bio, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_X_sim, 'c-', 'LineWidth', 2);
title('Biomass');
ylabel('c_X (g/L)');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w');
lgd = legend("Messung", "Simulation" );
lgd.TextColor = 'w';
grid on;

subplot(2, 1, 2);
plot(t_messung, y_glc, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_Glc_sim, 'y-', 'LineWidth', 2);
title('Glucose');
xlabel('BatchAge (h)');
ylabel('c_{Glc} (g/L)');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w');
lgd = legend("Messung", "Simulation" );
lgd.TextColor = 'w';
grid on;

set(gcf, 'Color', [0.2 0.2 0.2]);