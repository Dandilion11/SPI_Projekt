% Berechnung der Parameterunsicherheiten via Fisher-Informationsmatrix (FIM)
clear; clc; close all;

%% 1. Daten und optimierte Parameter laden
scriptDir = pwd;
loadPath_processed  = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');
load(loadPath_processed);

loadPath_processed  = fullfile(scriptDir, '..', 'Daten', 'p_opt.mat');
p_opt = (load(loadPath_processed));
p_opt = p_opt.p_opt;

projectRoot = pwd;
addpath(fullfile(projectRoot,'..','utils'),'-begin');


% Optimierte Parameter aus Aufgabe 4 (ValData = Experiment 04)
t_mes  = ValData.Biomasse.t;
y_bio  = ValData.Biomasse.y;
y_glc  = ValData.Glucose.y;
var_bio = ValData.Biomasse.var;
var_glc = ValData.Glucose.var;

x0 = [y_bio(1); y_glc(1)];

kinetic    = 3;
withOxygen = false;

nx = 2; np = 3;

%% 2. Kovarianzmatrix der Messung (invC)
% Messrauschen ist zeitabhaengig -> einmal gemittelt als konstante Matrix
C_mean = diag([mean(var_bio), mean(var_glc)]);
invC   = inv(C_mean);

%% 3. Erweiterte Simulation (Zustaende + Sensitivitaeten)
x0_ext = [x0; zeros(nx * np, 1)];  % XP(0) = 0

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[~, X_ext] = ode15s(@(t, x) Modell1_XP(t, x, p_opt, kinetic, withOxygen), ...
                     t_mes, x0_ext, options);

%% 4. FIM berechnen
FM = zeros(np, np);

% Messgleichung: h = x (Konzentrationen direkt gemessen)
dhdx = eye(nx);      % eye -> Einheitsmatrix der Größe nx
dhdp_direct = zeros(nx, np);

for k = 1:length(t_mes)
    x_k  = X_ext(k, 1:nx)';
    XP_k = reshape(X_ext(k, nx+1:end), nx, np);

    % dh/dp = dh/dx * XP  +  dh/dp_direkt
    dhdp = dhdx * XP_k + dhdp_direct;

    % Zeitpunktspezifische Messkovarianz
    C_k   = diag([var_bio(k), var_glc(k)]);
    invC_k = inv(C_k);

    % FIM aufaddieren
    FM = FM + dhdp' * invC_k * dhdp;
end

fprintf('Fisher-Informationsmatrix FM:\n');
disp(FM);

%% 5. Kovarianzmatrix und Parameteranalyse
CV = inv(FM);

[Corr, StdDev, relStdDev, EW, EV, CN] = Parameteranalyse(CV, p_opt);

fprintf('\n--- Parameterunsicherheiten ---\n');
fprintf('%-8s  %-10s  %-12s  %-14s\n', 'Param', 'Wert', 'StdAbw', 'Rel. Unsich.');
param_namen = {'mumax', 'KS', 'YXS'};
for i = 1:np
    fprintf('%-8s  %-10.4f  %-12.4f  %-14.1f%%\n', ...
        param_namen{i}, p_opt(i), StdDev(i), relStdDev(i)*100);
end
fprintf('\nKorrelationsmatrix:\n');
disp(Corr);
fprintf('Konditionszahl: %.1f\n', CN);

%% 6. Visualisierung der Parameterunsicherheit (Ellipsoid)
figure('Name', 'Parameterunsicherheit (FIM)', 'Position', [200, 200, 600, 500]);
SD = 1;  % 1 Standardabweichung ~ 65%
plot_gaussian_ellipsoid(p_opt, CV, SD);
xlabel('\mu_{max}');
ylabel('K_S');
zlabel('Y_{XS}');
title(sprintf('Parameterunsicherheit (1\\sigma)\nKonditionszahl: %.1f', CN));
grid on;