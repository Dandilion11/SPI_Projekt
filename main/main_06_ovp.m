%solves the optimal control problem to calculate optimal feeding trajectories (Task 6)
% main_06_ovp_modell3.m
clear; clc; close all;

%% 1. Laden der Ergebnisse aus der Parameteridentifikation (PI)
% Benötigt: p_opt_m3, x0_m3, invC (Messunsicherheit)
% Hier Platzhalter bis PI von Modell 3 vorliegt:
p_opt_m3 = zeros(12, 1); % Modell 3 hat 12 Basis-Parameter
x0_m3 = zeros(8, 1);     % Modell 3 hat 8 Zustände

%% 2. Konfiguration der Versuchsplanung (VP)
VP.t = 0:1:20;          % Zeitraster der Stellgrößenänderung und Probenahme (z.B. 20h)
VP.x0 = x0_m3;
VP.p = p_opt_m3;

% Initialisierung der u-Matrix für Modell 3 (10 Zeilen laut Modell3.m)
% Zeile 1: Zeitvektor
% Zeile 6: qGlc (Glucose-Feedrate) -> primäre Optimierungsvariable
% Zeile 7: cGlcF (Konzentration im Feed, ca. 450 g/L)
tu = VP.t;
u0_glc = 0.01 * ones(size(tu)); % Startschätzung für den Feed-Verlauf [L/h]

VP.u = zeros(10, length(tu));
VP.u(1, :) = tu;
VP.u(6, :) = u0_glc;
VP.u(7, :) = 450; % cGlcF konstant

%% 3. Grenzwertdefinitionen (Beschränkungen)
VP.CONS.umin = 0 * ones(size(u0_glc));
VP.CONS.umax = 0.05 * ones(size(u0_glc)); % Beispielgrenze für Pumpe

% Zustandsschranken [mX; mGlc; mNH4; mPO4; mEt; mB; V; DOT]
VP.CONS.xmin = [0; 0; 0; 0; 0; 0; 0.5; 0]; 
VP.CONS.xmax = [inf; inf; inf; inf; inf; inf; 2.0; inf]; % Reaktorvolumengrenze z.B. 2L

%% 4. Durchführung OVP
% Berechnung der initialen FIM des Vorversuchs (FM_old)
invC = eye(7); % Platzhalter für 7 Messgrößen aus Modell3_mgl.m
FM_old = zeros(length(p_opt_m3)); % Falls Vorversuchsdaten einfließen

fprintf('Starte Optimale Versuchsplanung für Modell 3...\n');
OVP_ERGEBNIS = OVP_Modell3(VP, FM_old, invC);

%% 5. Visualisierung des optimalen Feed-Profils
figure('Name', 'OVP Ergebnis Modell 3');
subplot(2,1,1)
stairs(OVP_ERGEBNIS.t, OVP_ERGEBNIS.u(6,:), 'LineWidth', 2);
ylabel('Optimaler Glucose-Feed q_{Glc} [L/h]');
grid on;

subplot(2,1,2)
plot(OVP_ERGEBNIS.t, OVP_ERGEBNIS.x(7,:), 'LineWidth', 2); % Reaktorvolumen
ylabel('Volumenverlauf V [L]');
xlabel('Zeit [h]');
grid on;