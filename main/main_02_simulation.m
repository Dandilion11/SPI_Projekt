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
    error('Data not found.\n');
else 
    fprintf('Data loaded:\n  %s\n',loadPath);
end 
load(loadPath);

t_bio   = TrainData.Biomasse.t;
y_bio   = TrainData.Biomasse.y;
var_bio = TrainData.Biomasse.var;

t_glc   = TrainData.Glucose.t;      
y_glc   = TrainData.Glucose.y;
var_glc = TrainData.Glucose.var;

t_o2 = TrainData.O2.t;
y_o2 = TrainData.O2.y;
var_o2 = TrainData.O2.var;

%% 2. Define Initial Conditions and Parameters
% Modell1 expects concentrations (cX, cGlc), not masses.
% Anfangswerte aus der jeweils ersten Messung.
x0 = [y_bio(1); y_glc(1)];

% Parameters for Monod: [mumax, KS, YXS] -> Startschätzwerte -> Richtige
% werden im nächstem Skript berechnet
p_guess = [0.3, 0.5, 0.15]; 

% Model Configuration Flags
kinetic = 3; % 3 = Monod
%--------------------------------------------------------------------------
withOxygen = true; % Exclude O2 state and parameters
%--------------------------------------------------------------------------
if withOxygen
    cO2_0    = 95;     % [%]     Start-DOT (aus O2-Messung zu Batch-Beginn)
    YXO      = 1.0;    % [gX/gO2] Ertragskoeffizient Biomasse/Sauerstoff (Schaetzwert)
    KLa      = 200;    % [1/h]    volumetrischer O2-Transferkoeffizient (Schaetzwert)
    cO2stern = 100;    % [%]      Saettigungs-DOT (Gleichgewicht mit Gasphase)

    x0      = [x0; cO2_0];                       % 3. Zustand ergaenzen
    p_guess = [p_guess, YXO, KLa, cO2stern];     % Parameter 4..6 ergaenzen
end



%% 3. Execute Simulation
% Simulation nur so lange laufen lassen wie Daten vorhanden sind:
t_end = max([t_bio(:)+1; t_glc(:)+1]);
t_sim = linspace(0, t_end, 200);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);


% Call Modell1 with the required flags
[~, X_sim] = ode45(@(t, x) Modell1(t, x, p_guess, kinetic, withOxygen), t_sim, x0, options);
fprintf('[OK] Simulation \n');

%% 4. Extract Data
% Division by V is no longer necessary as Modell1 outputs concentrations directly.
c_X_sim   = X_sim(:, 1);
c_Glc_sim = X_sim(:, 2);
if withOxygen
    cO2_sim = X_sim(:, 3);
end

%% 5. Plot Results
figure('Name', 'Simulation Test (Modell1)', 'Position', [200, 200, 900, 600]);
if withOxygen
    nRows = 3;
else 
    nRows = 2;
end
subplot(nRows,1,1)
% Messung mit Fehlerbalken -> +/- 1 Standardabweichung aus der Messvarianz
errorbar(t_bio, y_bio, sqrt(var_bio), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim, c_X_sim, 'LineWidth', 2);
title('Biomass');
ylabel('c_X (g/L)');
legend("Messung \pm \sigma", "Simulation",Location="southeast");
xlim([0 9.5])
grid on;

subplot(nRows, 1, 2);
%Glucose gegen eigenen Zeitvektor t_glc plotten + Standardabweichung aus
%Messung
errorbar(t_glc, y_glc, sqrt(var_glc), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim, c_Glc_sim, 'LineWidth', 2);
title('Glucose');
xlabel('BatchAge (h)');
ylabel('c_{Glc} (g/L)');
legend("Messung \pm \sigma", "Simulation",Location="northeast");
xlim([0 9.5])
grid on;

if withOxygen
    subplot(nRows, 1, 3);
    errorbar(t_o2, y_o2, sqrt(var_o2), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
    plot(t_sim, cO2_sim, 'LineWidth', 2); hold on;
    title('Sauerstoff (DOT)'); 
    ylabel('cO_2 (%)');
    
    legend("Messung \pm \sigma", "Simulation");
    xlim([0 9.5])
    grid on;
end
