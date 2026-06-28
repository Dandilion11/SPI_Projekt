% main_02_simulation_test.m
% Zweck: Vorwärts-Simulation des Modells mit geschätzten Parametern zur
% visuellen Überprüfung der Modellstruktur vor der Optimierung.

clear; clc; close all;

%% 1. Vorverarbeitete Daten laden
scriptDir = fileparts(mfilename('fullpath'));
loadPath = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');

% Überprüfen, ob Daten existieren
if ~isfile(loadPath)
    error('Daten nicht gefunden! Bitte zuerst main_01_data_preprocessing.m ausführen.');
end
load(loadPath);

% Trainingsdaten extrahieren (Experiment 03)
t_messung = TrainData.Biomasse.t;
y_bio     = TrainData.Biomasse.y;
y_glc     = TrainData.Glucose.y;
V         = TrainData.V; % Konstantes Reaktorvolumen/Gewicht für Batch-Phase

%% 2. Startbedingungen und Parameter definieren
% Anfangsmassen berechnen (Masse = Konzentration * Volumen)
m_X0   = y_bio(1) * V;
m_Glc0 = y_glc(1) * V;
x0     = [m_X0; m_Glc0];

% Erste Schätzung der biologischen Parameter für Hefe
% p = [mu_max (1/h), K_Glc (g/L), Y_Glc_X (g/g)]
p_guess = [0.3, 0.5, 0.15]; 

%% Modell-Auswahl: 1 = Monod, 2 = Tessier, 3 = Moser
model_type = 1;

%% Auswahl ob das Modell Sauerstoff mit einbeziehen soll oder nicht
withOxygen = false; 

%% 3. Simulation durchführen (Forward Problem)
% Zeitvektor für eine glatte Kurve (0 bis 20 Stunden, 200 Punkte)
t_sim = linspace(0, 20, 200);

% ODE Löser Optionen (Strenge Toleranzen für numerische Stabilität)
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

% Aufruf des ODE-Lösers
% Hinweis: ode_task2a_batch MUSS im Ordner /models/ liegen und der Pfad 
% muss in MATLAB bekannt sein (z.B. durch setup_project.m).
[~, X_sim] = ode15s(@(t, x) Modell1(t, x, p_guess, model_type, withOxygen), t_sim, x0, options);

%% 4. Post-Processing (Rückrechnung auf Konzentration)
% Die ODE rechnet mit absoluten Massen (g). Für den Plot müssen wir 
% diese zurück in Konzentrationen (g/L) wandeln (c = m/V).
c_X_sim   = X_sim(:, 1) / V;
c_Glc_sim = X_sim(:, 2) / V;

%% 5. Visuelle Überprüfung (Plot)
figure('Name', 'Simulation Test (Vorwärts-Simulation)', 'Position', [200, 200, 900, 600]);

% Biomasse
subplot(2, 1, 1);
plot(t_messung, y_bio, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_X_sim, 'c-', 'LineWidth', 2);
title(sprintf('Biomasse (Model %d) - mu_{max}=%.2f, K_S=%.2f, Y=%.2f', ...
    model_type, p_guess(1), p_guess(2), p_guess(3)), 'Color', 'w');
ylabel('c_X (g/L)', 'Color', 'w');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w');
legend('Messdaten (Exp 03)', 'Simulation (Geraten)', 'Location', 'NorthWest', 'TextColor', 'w');
grid on;

% Glucose
subplot(2, 1, 2);
plot(t_messung, y_glc, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6); hold on;
plot(t_sim, c_Glc_sim, 'y-', 'LineWidth', 2);
title('Glucose', 'Color', 'w');
xlabel('BatchAge (h)', 'Color', 'w');
ylabel('c_{Glc} (g/L)', 'Color', 'w');
set(gca, 'Color', [0.15 0.15 0.15], 'XColor', 'w', 'YColor', 'w');
legend('Messdaten (Exp 03)', 'Simulation (Geraten)', 'Location', 'NorthEast', 'TextColor', 'w');
grid on;

% Hintergrundfarbe der gesamten Figure anpassen
set(gcf, 'Color', [0.2 0.2 0.2]);