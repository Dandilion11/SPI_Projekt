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
t_bio  = ValData.Biomasse.t;
y_bio  = ValData.Biomasse.y;
var_bio = ValData.Biomasse.var;

t_glc  = ValData.Glucose.t;
y_glc  = ValData.Glucose.y;
var_glc = ValData.Glucose.var;

t_o2 = ValData.O2.t;
y_o2 = ValData.O2.y;  
var_o2 = ValData.O2.var;


kinetic    = 3; % kinetic 1 für Moser/Blachmann, Kinetic 2 für exponential und Kinetic 3
% Calculations only for Monod kinetic

withOxygen = true;

if withOxygen
    nx = 3; np = 5;
    x0 = [y_bio(1); y_glc(1); y_o2(1)];
    C_mean = diag([mean(var_bio), mean(var_glc), mean(var_o2)]); % Messrauschen ist zeitabhaengig -> einmal gemittelt als konstante Matrix
    t_all = unique([t_bio; t_glc; t_o2]);
else
    nx = 2; np = 3;
    x0 = [y_bio(1); y_glc(1)];
    C_mean = diag([mean(var_bio), mean(var_glc)]); % Messrauschen ist zeitabhaengig -> einmal gemittelt als konstante Matrix
    t_all = unique([t_bio; t_glc]);
end

%% 2. Kovarianzmatrix der Messung (invC)
invC   = inv(C_mean);

%% 3. Erweiterte Simulation (Zustaende + Sensitivitaeten)
x0_ext = [x0; zeros(nx * np, 1)];  % XP(0) = 0

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[~, X_ext] = ode15s(@(t, x) Modell1_XP_Monod(t, nx, np, x, p_opt, withOxygen), ...
                     t_all, x0_ext, options);

%% 4. FIM berechnen
FM = zeros(np, np);

% Sensitivitäten werden einzeln berechnet, weil o2 deutlich mehr
% Zeitschritte hat als die anderen Zustände (mehr Messungen)
% Am Ende werden die dann alle aufsummiert

[~, iBio] = ismember(t_bio, t_all);
[~, iGlc] = ismember(t_glc, t_all);

for k = 1:numel(iBio)
    XP_k = reshape(X_ext(iBio(k), nx+1:end), nx, np);

    s = XP_k(1, :).';              % Sensitivität Biomasse
    FM = FM + (s * s.') / var_bio(k);
end

for k = 1:numel(iGlc)
    XP_k = reshape(X_ext(iGlc(k), nx+1:end), nx, np);

    s = XP_k(2, :).';              % Sensitivität Glucose
    FM = FM + (s * s.') / var_glc(k);
end

if withOxygen
    [~, iO2] = ismember(t_o2, t_all);

    for k = 1:numel(iO2)
        XP_k = reshape(X_ext(iO2(k), nx+1:end), nx, np);

        s = XP_k(3, :).';          % Sensitivität Sauerstoff
        FM = FM + (s * s.') / var_o2(k);
    end
end

fprintf('Fisher-Informationsmatrix FM:\n');
disp(FM);
%% 5. Kovarianzmatrix und Parameteranalyse
CV = FM \ eye(np);

[Corr, StdDev, relStdDev, EW, EV, CN] = Parameteranalyse(CV, p_opt);

fprintf('\n--- Parameterunsicherheiten ---\n');
fprintf('%-8s  %-10s  %-12s  %-14s\n', 'Param', 'Wert', 'StdAbw', 'Rel. Unsich.');
if withOxygen
    param_namen = {'mumax', 'KS', 'YXS', 'YXO', 'KLa'};
else
    param_namen = {'mumax', 'KS', 'YXS'};
end
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
idxEll = [1 2 3];   % mumax, KS, YXS

plot_gaussian_ellipsoid(p_opt(idxEll), CV(idxEll, idxEll), SD);

xlabel('\mu_{max}');
ylabel('K_S');
zlabel('Y_{XS}');
title(sprintf('Parameterunsicherheit (1\\sigma)\nKonditionszahl: %.1f', CN));
grid on;