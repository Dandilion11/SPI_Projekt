% main_02_simulation_test.m
clear; clc; close all;

% 1. Load Preprocessed Data
scriptDir = pwd;
loadPath = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_FedBatch_Modell3.mat');
projectRoot = pwd;
addpath(fullfile(projectRoot, '..','Modelle'),'-begin');

rehash;
clear Modell3

if ~isfile(loadPath)
    error('Data not found.\n');
else 
    fprintf('Data loaded:\n  %s\n',loadPath);
end 
load(loadPath);

Data = TrainData; % Geht auch ValData -> zum schnelleren wechseln

% Anfangswerte und Parameter
%   [mumax,KS,YXS,YAmX,YPhX,YB_Am, KLa, YXO]
p_guess = [0.3, 0.5, 0.15, 0.05, 0.02, 1.0, 200, 1.0];
x0 = Data.x0;
u = Data.u;

% Simulation
t_start = Data.u(1,1);              % = 0.03
t_end   = max(Data.Biomasse.t(:)) + 1;
t_sim   = linspace(t_start, t_end, 300);

fprintf('u-Startzeit = %.4f, t_sim-Startzeit = %.4f\n', Data.u(1,1), t_sim(1));
assert(size(Data.u,1) >= 10, 'u-Matrix hat zu wenige Zeilen (%d)', size(Data.u,1));

DOTstern = max(Data.O2.y);

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[~, X3]  = ode45(@(t, x) Modell3_Jannis(t, x, u, p_guess, DOTstern), t_sim, x0, options);

V    = X3(:,1);
cX   = X3(:,2)./V;   % Biomasse [g/L]
cGlc = X3(:,3)./V;   % Glucose  [g/L]
cAm  = X3(:,4)./V;   % Ammonium [g/L]
cPh  = X3(:,5)./V;   % Phosphat [g/L]
mB   = X3(:,6);      % kumulierte Base
DOT  = X3(:,7);      % Sauerstoff [%]

% Plot Modell3
figure('Name','Simulation Modell3 (Fed-Batch + Ethanol)','Position',[200 40 950 1000]);

subplot(7,1,1);
errorbar(Data.Biomasse.t, Data.Biomasse.y, sqrt(Data.Biomasse.var), 'o', ...
         'MarkerFaceColor','b','MarkerSize',4); hold on;
plot(t_sim, cX, 'LineWidth', 2);
title('Modell3 - Biomasse'); ylabel('c_X (g/L)');
legend("Messung \pm \sigma","Simulation"); grid on;

subplot(7,1,2);
errorbar(Data.Glucose.t, Data.Glucose.y, sqrt(Data.Glucose.var), 'o', ...
         'MarkerFaceColor','b','MarkerSize',4); hold on;
plot(t_sim, cGlc, 'LineWidth', 2);
title('Modell3 - Glucose'); ylabel('c_{Glc} (g/L)');
legend("Messung \pm \sigma","Simulation"); grid on;

subplot(7,1,3);
errorbar(Data.Ammonium.t, Data.Ammonium.y, sqrt(Data.Ammonium.var), 'o', ...
         'MarkerFaceColor','b','MarkerSize',4); hold on;
plot(t_sim, cAm, 'LineWidth', 2);
title('Modell3 - Ammonium'); ylabel('c_{Am} (g/L)');
legend("Messung \pm \sigma","Simulation"); grid on;

subplot(7,1,4);
errorbar(Data.Phosphat.t, Data.Phosphat.y, sqrt(Data.Phosphat.var), 'o', ...
         'MarkerFaceColor','b','MarkerSize',4); hold on;
plot(t_sim, cPh, 'LineWidth', 2);
title('Modell3 - Phosphat'); ylabel('c_{Ph} (g/L)');
legend("Messung \pm \sigma","Simulation"); grid on;

subplot(7,1,5);
errorbar(Data.Base.t, Data.Base.y, sqrt(Data.Base.var), 'o', ...
         'MarkerFaceColor','b','MarkerSize',4); hold on;
plot(t_sim, mB, 'LineWidth', 2);
title('Modell3 - Base'); ylabel('m_B');
legend("Messung \pm \sigma","Simulation"); grid on;

subplot(7,1,6);
errorbar(Data.O2.t, Data.O2.y, sqrt(Data.O2.var), 'o', ...
         'MarkerFaceColor','b','MarkerSize',4); hold on;
plot(t_sim, DOT, 'LineWidth', 2);
title('Modell3 - Sauerstoff (DOT)'); ylabel('DOT (%)');
legend("Messung \pm \sigma","Simulation"); grid on;






