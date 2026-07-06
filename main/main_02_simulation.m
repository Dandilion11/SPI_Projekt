% main_02_simulation_test.m
clear; clc; close all;

%% 1. Load Preprocessed Data
scriptDir = pwd;
loadPath = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');
projectRoot = pwd;
addpath(fullfile(projectRoot, '..','Modelle'),'-begin');

rehash;
clear Modell1 Modell2

if ~isfile(loadPath)
    error('Data not found.\n');
else 
    fprintf('Data loaded:\n  %s\n',loadPath);
end 
load(loadPath);

Data = TrainData; % Geht auch ValData -> zum schnelleren wechseln

t_bio   = Data.Biomasse.t;
y_bio   = Data.Biomasse.y;
var_bio = Data.Biomasse.var;

t_glc   = Data.Glucose.t;      
y_glc   = Data.Glucose.y;
var_glc = Data.Glucose.var;

t_o2 = Data.O2.t;
y_o2 = Data.O2.y;
var_o2 = Data.O2.var;

%% A2. Anfangswerte und Parameter
kinetic = 3; % 3 = Monod
%--------------------------------------------------------------------------
withOxygen = true; 
%--------------------------------------------------------------------------
% Anfangswerte aus der jeweils ersten Messung.
x0 = [y_bio(1); y_glc(1)];

% Parameters for Monod: [mumax, KS, YXS] -> Startschätzwerte -> Richtige
% werden im nächstem Skript berechnet
p_guess = [0.3, 0.5, 0.15]; 

% Model Configuration Flags

if withOxygen
    cO2_0    = y_o2(1);     % [%]     Start-DOT (aus O2-Messung zu Batch-Beginn)
    YXO      = 1.0;    % [gX/gO2] Ertragskoeffizient Biomasse/Sauerstoff (Schaetzwert)
    KLa      = 200;    % [1/h]    volumetrischer O2-Transferkoeffizient (Schaetzwert)
    cO2stern = 100;    % [%]      Saettigungs-DOT (Gleichgewicht mit Gasphase)

    x0      = [x0; cO2_0];                       % 3. Zustand ergaenzen
    p_guess = [p_guess, YXO, KLa, cO2stern];     % Parameter 4..6 ergaenzen
end



%% A3. Execute Simulation
% Simulation nur so lange laufen lassen wie Daten vorhanden sind:
t_end = max([t_bio(:)+1; t_glc(:)+1]);
t_sim = linspace(0, t_end, 200);
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);


% Call Modell1 with the required flags
[~, X_sim] = ode45(@(t, x) Modell1(t, x, p_guess, kinetic, withOxygen), t_sim, x0, options);
fprintf('[OK] Simulation \n');

%% A4. Extract Data
% Division by V is no longer necessary as Modell1 outputs concentrations directly.
c_X_sim   = X_sim(:, 1);
c_Glc_sim = X_sim(:, 2);
if withOxygen
    cO2_sim = X_sim(:, 3);
end

%% A5. Plot Results
figure('Name', 'Simulation Test (Modell1)', 'Position', [200, 200, 900, 600]);
if withOxygen
    nRows = 3;
else 
    nRows = 2;
end
subplot(nRows,1,1)
% Messung mit Fehlerbalken -> +/- 1 Standardabweichung aus der Messvarianz
errorbar(t_bio, y_bio, sqrt(var_bio), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim, c_X_sim, 'LineWidth', 2);
title('Biomass');
ylabel('c_X (g/L)');
legend("Messung \pm \sigma", "Simulation",Location="southeast");
xlim([0 9.5])
grid on;

subplot(nRows, 1, 2);
%Glucose gegen eigenen Zeitvektor t_glc plotten + Standardabweichung aus
%Messung
errorbar(t_glc, y_glc, sqrt(var_glc), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim, c_Glc_sim, 'LineWidth', 2);
title('Glucose');
xlabel('BatchAge (h)');
ylabel('c_{Glc} (g/L)');
legend("Messung \pm \sigma", "Simulation",Location="northeast");
xlim([0 9.5])
grid on;

if withOxygen
    subplot(nRows, 1, 3);
    errorbar(t_o2, y_o2, sqrt(var_o2), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
    plot(t_sim, cO2_sim, 'LineWidth', 2); hold on;
    title('Sauerstoff (DOT)'); 
    ylabel('cO_2 (%)');
    
    legend("Messung \pm \sigma", "Simulation");
    xlim([0 9.5])
    grid on;
end




%% MODELL 2 Geht ab hier los 
%% B6. Anfangswerte und Parameter
clear Modell2

x0_m2 = [Data.Biomasse.y(1); Data.Glucose.y(1); Data.Ammonium.y(1); Data.Base.y(1); Data.O2.y(1)];

% Startschaetzwerte: GEGENCHECKEN!!!!!
%   YBam = Yield Base/Ammonium, YAmX = Yield Ammonium/Biomasse
p_guess_m2 = [0.3, 0.5, 0.15, 1.0, 0.05, 1.0, 200, 100];

%% B7. Simulation
t_end_m2 = max([Data.Biomasse.t(:)+1; Data.Glucose.t(:)+1; Data.Ammonium.t(:)+1; Data.Base.t(:)+1; Data.O2.t(:)+1]);

t_sim_m2 = linspace(0, t_end_m2, 300);

[~, X2]  = ode45(@(t, x) Modell2(t, x, p_guess_m2, kinetic), t_sim_m2, x0_m2, options);


%% B8. Plot Modell2
figure('Name','Simulation Test (Modell2: +Base +O2)','Position',[250 80 950 950]);

subplot(5,1,1);
errorbar(Data.Biomasse.t, Data.Biomasse.y, sqrt(Data.Biomasse.var), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:,1), 'LineWidth', 2);
title('Modell2 - Biomasse'); 
ylabel('c_X (g/L)'); 
legend("Messung \pm \sigma","Simulation"); 
xlim([0 10])
grid on;

subplot(5,1,2);
errorbar(Data.Glucose.t, Data.Glucose.y, sqrt(Data.Glucose.var), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:,2), 'LineWidth', 2);
title('Modell2 - Glucose'); 
ylabel('c_{Glc} (g/L)');  
xlim([0 10])
legend("Messung \pm \sigma","Simulation");
grid on;

subplot(5,1,3);
errorbar(Data.Ammonium.t, Data.Ammonium.y, sqrt(Data.Ammonium.var), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:,3), 'LineWidth', 2);
title('Modell2 - Ammonium'); 
ylabel('c_{Am} (g/L)'); 
xlim([0 10])
legend("Messung \pm \sigma","Simulation");
grid on;

subplot(5,1,4);
errorbar(Data.Base.t, Data.Base.y, sqrt(Data.Base.var), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:,4), 'LineWidth', 2);
title('Modell2 - Base'); 
ylabel('c_{Base}'); 
xlim([0 10])
legend("Messung \pm \sigma","Simulation");
grid on;

subplot(5,1,5);
errorbar(Data.O2.t, Data.O2.y, sqrt(Data.O2.var), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:,5), 'LineWidth', 2);
title('Modell2 - Sauerstoff (DOT)'); 
ylabel('cO_2 (%)'); 
xlabel('BatchAge (h)');
xlim([0 10])
legend("Messung \pm \sigma","Simulation");
grid on;






