%solves the optimal control problem to calculate optimal feeding trajectories (Task 6)
% main_06_ovp_modell3.m
clear; clc; close all;

%% 1. Laden der Ergebnisse aus der Parameteridentifikation (PI)
% Benötigt: p_opt_m3, x0_m3, invC (Messunsicherheit)
% Hier Platzhalter bis PI von Modell 3 vorliegt:
load(fullfile(pwd,'Daten','p_opt', 'p_opt_Modell3_woEtOH.mat'));
load(fullfile(pwd,'..', 'Daten','Daten_Processed','Processed_FedBatch_Modell3.mat'))
addpath (fullfile(pwd,'..','\utils'))
p_opt_m3 = p_opt;
x0_m3    = TrainData.x0;  

%% 2. Konfiguration der Versuchsplanung (VP)
ab = [0.02 0.015; 0.06 0.25; 0.06 0.01; 0.07 0.01; 0.01 0.01; 0.02 0.5];
c_nom = [5; 2; 0.5; 0.2; 0.05; 60]; % repräsentative Konz. [g/L, g/L, g/L, g/L, L, %]
sigma_nom = (ab(:,1).*c_nom + ab(:,2)).^2;
C    = diag(sigma_nom);
invC = inv(C);
VP.t = 0:20; % 
tu   = VP.t;                  % Stellgroessenraster identisch
VP.x0 = x0_m3;
VP.p  = p_opt_m3;
VP.DOTstern = max(TrainData.O2.y);

u0_glc  = 0.05 * ones(size(tu));   % realistischerer Startwert
u_am_fix = 0.01 * ones(size(tu));  % proportionaler Am-Feed
u_ph_fix = 0.002 * ones(size(tu)); % proportionaler Ph-Feed

VP.u = zeros(10, length(tu));
VP.u(1, :) = tu;
VP.u(2, :) = u_am_fix;  VP.u(3, :) = 80;   % Am-Feed [L/h], [g/L]
VP.u(4, :) = u_ph_fix;  VP.u(5, :) = 40;   % Ph-Feed [L/h], [g/L]
VP.u(6, :) = u0_glc;    VP.u(7, :) = 450;  % Glc-Feed [L/h], [g/L]

%% 3. Grenzwertdefinitionen (Beschränkungen)
VP.CONS.umin = zeros(size(u0_glc));
VP.CONS.umax = 0.20 * ones(size(u0_glc)); % max Feed [L/h]

% Zustandsschranken [V; mX; mGlc; mAm; mPh; mB; DOT]
VP.CONS.xmin = [0.5; 0; 0; 0; 0; 0; 0];
VP.CONS.xmax = [14.0; inf; inf; inf; inf; inf; 100];


%% 4. Durchführung OVP
% Berechnung der initialen FIM des Vorversuchs (FM_old)
[FM_old, ~, ~] = simulation_fim_modell3(VP.t, VP.x0, VP.u, VP.p, invC, VP.DOTstern);

FM_norm_alt = diag(VP.p) * FM_old * diag(VP.p);
CV_alt      = inv(FM_norm_alt);
[Corr, StdDev, relStdDev, ~, ~, CN] = Parameteranalyse(CV_alt, VP.p);
fprintf('Konditionszahl vor OVP: %.1f\n', CN);

figure; plot_gaussian_ellipsoid(VP.p(1:3), CV_alt(1:3,1:3), 1);
title('Parameterunsicherheit vor OVP');

fprintf('Starte Optimale Versuchsplanung für Modell 3...\n');
OVP_ERGEBNIS = OVP_Modell3(VP, FM_old, invC);

[FM_neu, ~, ~] = simulation_fim_modell3(VP.t, VP.x0, ...
                     OVP_ERGEBNIS.u, VP.p, invC, VP.DOTstern);

FM_total      = FM_old + FM_neu;
FM_norm_neu   = diag(VP.p) * FM_total * diag(VP.p);
CV_neu        = inv(FM_norm_neu);
[~, StdDev2, relStdDev2, ~, ~, CN2] = Parameteranalyse(CV_neu, VP.p);

figure; plot_gaussian_ellipsoid(VP.p(1:3), CV_neu(1:3,1:3), 1);
title('Erwartete Parameterunsicherheit nach OVP');
%% 5. Visualisierung des optimalen Feed-Profils
figure('Name', 'OVP Ergebnis Modell 3');
subplot(2,1,1)
stairs(OVP_ERGEBNIS.t, OVP_ERGEBNIS.u(6,:), 'LineWidth', 2);
ylabel('Optimaler Glucose-Feed q_{Glc} [L/h]');
grid on;

subplot(2,1,2)
plot(OVP_ERGEBNIS.t, OVP_ERGEBNIS.x(1,:), 'LineWidth', 2); % Reaktorvolumen
ylabel('Volumenverlauf V [L]');
xlabel('Zeit [h]');
grid on;
