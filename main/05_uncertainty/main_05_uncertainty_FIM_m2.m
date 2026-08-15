% main_05_uncertainty_FIM_m2.m
% Berechnung der Parameterunsicherheiten via FIM für MODELL 2 (+Base, +O2)
clear; clc; close all;

%% 1. Daten laden
scriptDir = pwd;
loadPath_processed  = fullfile(scriptDir,'..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');
load(loadPath_processed);

projectRoot = pwd;
addpath(fullfile(projectRoot,'..','utils'),'-begin');
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');

% Validierungsdaten extrahieren
Data = TrainData;
t_bio = Data.Biomasse.t(:);   y_bio = Data.Biomasse.y(:);   var_bio = Data.Biomasse.var(:);
t_glc = Data.Glucose.t(:);    y_glc = Data.Glucose.y(:);    var_glc = Data.Glucose.var(:);
t_am  = Data.Ammonium.t(:);   y_am  = Data.Ammonium.y(:);   var_am  = Data.Ammonium.var(:);
t_ba  = Data.Base.t(:);       y_ba  = Data.Base.y(:);       var_ba  = Data.Base.var(:);
t_o2  = Data.O2.t(:);         y_o2  = Data.O2.y(:);         var_o2  = Data.O2.var(:);

t_mes = unique([t_bio; t_glc; t_am; t_ba; t_o2]);

% Startwerte für 5 Zustände
x0 = [y_bio(1); y_glc(1); y_am(1); y_ba(1); y_o2(1)];

% Optimierte Parameter aus PI
% Reihenfolge: [mumax, KS, YXS, YBam, YAmX, YXO, KLa]
loadPath_processed  = fullfile(scriptDir, '..', 'Daten', 'p_opt_Modell2.mat');
p_opt = (load(loadPath_processed));
p_opt_m2 = p_opt.p_opt_m2;

kinetic = 3; % Monod
nx = 5; 
np = 7;

%% 2. Kovarianzmatrix der Messung (invC)
% Messrauschen gemittelt als konstante 5x5 Diagonalmatrix
C_mean = diag([mean(var_bio), mean(var_glc), mean(var_am), mean(var_ba), mean(var_o2)]);
invC   = inv(C_mean);
DOTstern = max(y_o2);

%% 3. Erweiterte Simulation (Zustände + Sensitivitäten)
x0_ext = [x0; zeros(nx * np, 1)];  

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[~, X_ext] = ode15s(@(t, x) Modell2_XP_Monod(t, x, p_opt_m2, DOTstern), t_mes, x0_ext, options);

%% 4. FIM berechnen
FM = zeros(np, np);

t_list   = {t_bio, t_glc, t_am, t_ba, t_o2};
var_list = {var_bio, var_glc, var_am, var_ba, var_o2};

% Zustandsindex der jeweiligen Messgröße:
% Biomasse -> x(1), Glucose -> x(2), Ammonium -> x(3), Base -> x(4), O2 -> x(5)
state_idx = [1, 2, 3, 4, 5];

for m = 1:nx
    t_m   = t_list{m};
    var_m = var_list{m};

    [tf, idx_mes] = ismember(t_m, t_mes);

    if any(~tf)
        error('Nicht alle Messzeitpunkte wurden in t_mes gefunden.');
    end

    for j = 1:numel(t_m)
        k = idx_mes(j);

        XP_k = reshape(X_ext(k, nx+1:end), nx, np);

        % Sensitivität der gemessenen Größe nach Parametern
        s = XP_k(state_idx(m), :).';   % np x 1

        % FIM-Beitrag: s*s' / Varianz
        FM = FM + (s * s.') / var_m(j);
    end
end

fprintf('Fisher-Informationsmatrix FM (Modell 2):\n');
disp(FM);
fprintf('Rang(FM) = %d von %d\n', rank(FM), np);
fprintf('rcond(FM) = %.3e\n', rcond(FM));

%% 5. Kovarianzmatrix und Parameteranalyse
CV = FM \ eye(np);

[Corr, StdDev, relStdDev, EW, EV, CN] = Parameteranalyse(CV, p_opt_m2);

fprintf('\n--- Parameterunsicherheiten (Modell 2) ---\n');
fprintf('%-8s  %-10s  %-12s  %-14s\n', 'Param', 'Wert', 'StdAbw', 'Rel. Unsich.');
param_namen = {'mumax', 'KS', 'YXS', 'YBam', 'YAmX', 'YXO', 'KLa'};
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