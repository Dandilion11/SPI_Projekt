% main_05_uncertainty_FIM_m2.m
% Berechnung der Parameterunsicherheiten via FIM für MODELL 2 (+Base, +O2)
clear; clc; close all;

%% 1. Daten laden
scriptDir = pwd;
loadPath_processed  = fullfile(scriptDir, 'Matlab_Code', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');
load(loadPath_processed);

projectRoot = pwd;
addpath(fullfile(projectRoot,'..','utils'),'-begin');
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');

% Validierungsdaten extrahieren
Data = ValData;
t_mes  = Data.Biomasse.t;
y_bio  = Data.Biomasse.y;   var_bio = Data.Biomasse.var;
y_glc  = Data.Glucose.y;    var_glc = Data.Glucose.var;
y_am   = Data.Ammonium.y;   var_am  = Data.Ammonium.var;
y_ba   = Data.Base.y;       var_ba  = Data.Base.var;
y_o2   = Data.O2.y;         var_o2  = Data.O2.var;

% Startwerte für 5 Zustände
x0 = [y_bio(1); y_glc(1); y_am(1); y_ba(1); y_o2(1)];

% Optimierte Parameter aus PI
% Reihenfolge: [mumax, KS, YXS, YBam, YAmX, YXO, KLa, cO2stern]
p_opt_m2 = [0.3795; 5.0000; 0.1473; 1.2283; 0.0554; 2.6039; 288.1887; 76.8628]; 

kinetic = 3; % Monod
nx = 5; 
np = 8;

%% 2. Kovarianzmatrix der Messung (invC)
% Messrauschen gemittelt als konstante 5x5 Diagonalmatrix
C_mean = diag([mean(var_bio), mean(var_glc), mean(var_am), mean(var_ba), mean(var_o2)]);
invC   = inv(C_mean);

%% 3. Erweiterte Simulation (Zustände + Sensitivitäten)
x0_ext = [x0; zeros(nx * np, 1)];  

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[~, X_ext] = ode15s(@(t, x) Modell2_XP(t, x, p_opt_m2, kinetic), t_mes, x0_ext, options);

%% 4. FIM berechnen
FM = zeros(np, np);

dhdx = eye(nx);      
dhdp_direct = zeros(nx, np);

for k = 1:length(t_mes)
    XP_k = reshape(X_ext(k, nx+1:end), nx, np);
    dhdp = dhdx * XP_k + dhdp_direct;

    C_k   = diag([var_bio(k), var_glc(k), var_am(k), var_ba(k), var_o2(k)]);
    invC_k = inv(C_k);

    FM = FM + dhdp' * invC_k * dhdp;
end

fprintf('Fisher-Informationsmatrix FM (Modell 2):\n');
disp(FM);

%% 5. Kovarianzmatrix und Parameteranalyse
CV = inv(FM);

[Corr, StdDev, relStdDev, EW, EV, CN] = Parameteranalyse(CV, p_opt_m2);

fprintf('\n--- Parameterunsicherheiten (Modell 2) ---\n');
fprintf('%-8s  %-10s  %-12s  %-14s\n', 'Param', 'Wert', 'StdAbw', 'Rel. Unsich.');
param_namen = {'mumax', 'KS', 'YXS', 'YBam', 'YAmX', 'YXO', 'KLa', 'cO2stern'};
for i = 1:np
    fprintf('%-8s  %-10.4f  %-12.4f  %-14.1f%%\n', ...
        param_namen{i}, p_opt_m2(i), StdDev(i), relStdDev(i)*100);
end

fprintf('\nKonditionszahl: %.1f\n', CN);

%% 6. Visualisierung der Parameterunsicherheit (3D-Ellipsoid)
figure('Name', 'Parameterunsicherheit Modell 2 (FIM)', 'Position', [200, 200, 600, 500]);
SD = 1;  % 1 Standardabweichung ~ 65%

p_plot = p_opt_m2(1:3);
CV_plot = CV(1:3, 1:3);

plot_gaussian_ellipsoid(p_plot, CV_plot, SD);
xlabel('\mu_{max}');
ylabel('K_S');
zlabel('Y_{XS}');
title(sprintf('Parameterunsicherheit (1\\sigma) Modell 2\nKonditionszahl: %.1f', CN));
grid on;