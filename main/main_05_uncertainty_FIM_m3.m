clear; clc; close all;

projectRoot = pwd;
load(fullfile(projectRoot,'..','Daten','Daten_Processed','Processed_FedBatch_Modell3.mat'));
load(fullfile(projectRoot,'..','Daten','p_opt_Modell3_woEtOH.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');

Data  = TrainData;
u     = Data.u;
p_opt_m3 = p_opt;   % [mumax, KS, YXS, YAmX, YPhX, YB_Am, KLa, YXO]
param_namen = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO'};
np = 8;  nx = 7;  ny = 6;

DOTstern = 100;  % Enfors Gl. 6.15 / Krämer 2017

% --- invC: zustandsabhängige Messkovarianz (Krämer 2016, Tab. 2) ---
% Repräsentative Konzentrationen aus Trainingsdaten
c_rep = [mean(Data.Biomasse.y), mean(Data.Glucose.y), ...
         mean(Data.Ammonium.y), mean(Data.Phosphat.y), ...
         mean(Data.Base.y),     mean(Data.O2.y)];
ab = [0.02 0.015; 0.06 0.25; 0.06 0.01; 0.07 0.01; 0.01 10; 0.02 0.5];
sigma_rep = (ab(:,1).*c_rep(:) + ab(:,2)).^2;
C    = diag(sigma_rep);
invC = inv(C);

% --- FIM und Simulation ---
t_fim = unique([Data.Biomasse.t; Data.Glucose.t; Data.Ammonium.t; ...
                Data.Phosphat.t; Data.Base.t;    Data.O2.t]);
[FM, ~, ~] = simulation_fim_modell3(t_fim, Data.x0, u, p_opt_m3, invC, DOTstern);

% --- Normierte FIM (SPI-Skript Gl. 5.17) ---
FM_norm = diag(p_opt_m3) * FM * diag(p_opt_m3);
CV      = inv(FM_norm);

[Corr, StdDev, relStdDev, EW, EV, CN] = Parameteranalyse(CV, p_opt_m3);

fprintf('\n--- FIM Parameterunsicherheiten (Modell3_woEtOH) ---\n');
fprintf('%-10s  %-10s  %-12s  %-14s\n','Param','Wert','StdAbw','Rel.[%%]');
for i = 1:np
    fprintf('%-10s  %-10.4f  %-12.4f  %-14.1f\n', ...
        param_namen{i}, p_opt_m3(i), StdDev(i), relStdDev(i)*100);
end
fprintf('Konditionszahl: %.1f\n', CN);

% --- 3D-Ellipsoid ---
figure('Name','Parameterunsicherheit Modell3_woEtOH (FIM)');
plot_gaussian_ellipsoid(p_opt_m3(1:3), CV(1:3,1:3), 1);
xlabel('\mu_{max}'); ylabel('K_S'); zlabel('Y_{XS}');
title(sprintf('1\\sigma-Ellipsoid Modell3\\_woEtOH  (CN = %.1f)', CN));
grid on;