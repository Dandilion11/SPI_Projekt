clear; clc; close all;

% Datei laden (Beispiel)
scriptDir = fileparts(mfilename('fullpath'));
pfad = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef05.mat');
daten = load(pfad);
u_matrix = daten.Mess.u;

% 1. Zeitvektor extrahieren (Zeile 1)
t_u = u_matrix(1, :);

% 2. Volumenströme extrahieren (Einheit: L/h)
u_Am   = u_matrix(2, :); % Ammonium-Feed
u_Ph   = u_matrix(4, :); % Phosphat-Feed
u_Glc  = u_matrix(6, :); % Glucose-Feed
u_Acid = u_matrix(9, :); % Säure
u_Base = u_matrix(10,:); % Base

% Konstante Konzentrationen (Zur Verifizierung, nicht geplottet)
c_Am_in  = u_matrix(3, 1); % 30 g/L
c_Ph_in  = u_matrix(5, 1); % 24 g/L
c_Glc_in = u_matrix(7, 1); % 450 g/L

% 3. Plotten
figure('Name', 'Stellgrößen (Volumenströme)', 'Position', [100, 100, 800, 600]);

subplot(3,2,1);
stairs(t_u, u_Am, 'LineWidth', 1.5);
title('Ammonium Feed (u_{Am})');
ylabel('L/h'); grid on;

subplot(3,2,2);
stairs(t_u, u_Ph, 'LineWidth', 1.5);
title('Phosphat Feed (u_{Ph})');
ylabel('L/h'); grid on;

subplot(3,2,3);
stairs(t_u, u_Glc, 'LineWidth', 1.5);
title('Glucose Feed (u_{Glc})');
ylabel('L/h'); grid on;

subplot(3,2,4);
stairs(t_u, u_Acid, 'LineWidth', 1.5);
title('Acid Feed (u_{Acid})');
ylabel('L/h'); grid on;

subplot(3,2,5);
stairs(t_u, u_Base, 'LineWidth', 1.5);
title('Base Feed (u_{Base})');
xlabel('BatchAge (h)'); ylabel('L/h'); grid on;

sgtitle(sprintf('Eingangsvektor u - %s', daten.Mess.Name));