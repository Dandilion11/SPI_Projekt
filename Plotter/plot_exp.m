% main_00_data_visualization.m
clear; clc; close all;

%% 1. Configuration
experiment_nummer = 3; % Change this integer to load different experiments (e.g., 3, 4)

% Error model parameters (Krämer & King 2017)
a_bio = 0.02; b_bio = 0.015;
a_glc = 0.06; b_glc = 0.25;

%% 2. Data Loading
%scriptDir = fileparts(mfilename('fullpath'));
scriptDir = pwd;
dateiname = sprintf('Mess_RamScDef%02d.mat', experiment_nummer);
pfad = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', dateiname);

if ~isfile(pfad)
    error('File %s not found. Verify the directory structure.', dateiname);
end

daten = load(pfad);
Mess = daten.Mess;

%% 3. Data Extraction & Error Calculation
% Offline Data (At-line data like BiomasseOD/GlucoseBrix are explicitly ignored)
t_bio = Mess.Messdaten.Biomasse.BatchAge;
y_bio = Mess.Messdaten.Biomasse.Wert;
sigma_bio = a_bio .* y_bio + b_bio; % Standard deviation

t_glc = Mess.Messdaten.Glucose.BatchAge;
y_glc = Mess.Messdaten.Glucose.Wert;
sigma_glc = a_glc .* y_glc + b_glc;

% Input Data (u Matrix - L/h)
t_u    = Mess.u(1, :);
u_Ph  = Mess.u(2, :);  
u_Am   = Mess.u(4, :);  % Ammonium is row 4
u_Glc  = Mess.u(6, :);  
u_Acid = Mess.u(9, :);  % Assuming Acid remains row 9
u_Base = Mess.u(10, :); % Base is row 10

%% 4. Plotting (Updated for Dark Themes)
figure('Name', sprintf('Experiment %02d Data', experiment_nummer), 'Position', [100, 100, 1000, 800]);

% Subplot 1: Biomass
subplot(3, 2, 1);
% Changed 'ko' (black) to 'wo' (white) and MarkerFaceColor to 'w'
errorbar(t_bio, y_bio, sigma_bio, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.2);
xline(20, 'r--', 'Batch Limit (Task 2a)', 'LabelVerticalAlignment', 'bottom');
title('Biomass Concentration (Raw)');
xlabel('BatchAge (h)');
ylabel('c_X (g/L)');
grid on;

% Subplot 2: Glucose
subplot(3, 2, 2);
errorbar(t_glc, y_glc, sigma_glc, 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 6, 'LineWidth', 1.2);
xline(20, 'r--', 'Batch Limit (Task 2a)', 'LabelVerticalAlignment', 'bottom');
title('Glucose Concentration (Raw)');
xlabel('BatchAge (h)');
ylabel('c_{Glc} (g/L)');
grid on;

% Subplot 3: Volumetric Feeds (Nutrients)
subplot(3, 2, [3, 4]);
stairs(t_u, u_Glc, 'b', 'LineWidth', 1.5); hold on;
stairs(t_u, u_Am, 'g', 'LineWidth', 1.5);
stairs(t_u, u_Ph, 'm', 'LineWidth', 1.5);
title('Volumetric Feed Rates (Inputs)');
xlabel('BatchAge (h)');
ylabel('Flow Rate (L/h)');
legend('Glucose (u_{Glc})', 'Ammonium (u_{Am})', 'Phosphate (u_{Ph})', 'Location', 'NorthWest');
grid on;

% Subplot 4: Volumetric Feeds (pH Control)
subplot(3, 2, [5, 6]);
% Changed 'k' (black) to 'y' (yellow) for Base feed
stairs(t_u, u_Base, 'y', 'LineWidth', 1.5); hold on;
stairs(t_u, u_Acid, 'r', 'LineWidth', 1.5);
title('pH Control Volumes');
xlabel('BatchAge (h)');
ylabel('Flow Rate (L/h)');
legend('Base (u_{Base})', 'Acid (u_{Acid})', 'Location', 'NorthWest', 'TextColor', 'w');
grid on;

sgtitle(sprintf('Dataset Visualization: RamScDef%02d', experiment_nummer), 'Color', 'w');