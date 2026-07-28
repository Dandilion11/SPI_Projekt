% main_04_parameter_fitting.m
clear; clc; close all;
%{
%% 1. Load Preprocessed Data
load('Processed_Batch_Data.mat');
t_messung = TrainData.Biomasse.t;
y_bio = TrainData.Biomasse.y;
y_glc = TrainData.Glucose.y;
var_bio = TrainData.Biomasse.var;
var_glc = TrainData.Glucose.var;
%}

% Exp3 = TrainData
% Exp4 = ValData
% mit ValData läuft die PI
% Ich glaube bei Exp3 wird Glucose zugeführt (ca. stunde 8). Dadurch wäre
% das kein geeigenetes experiment mehr für das simpe model.
% Wenn tatsächlich glucose gefüttert wird, könnten wir bei exp3 die daten
% bis zum füttern nehmen

%% --------------------------------------------------------------------------
%% MODELL 1 -> Visualisierung selber aufrufen
%% --------------------------------------------------------------------------
% 1. Load Preprocessed Data
projectRoot = pwd;
load(fullfile(projectRoot,'..','Daten/Daten_Processed/Processed_Batch_Data.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');

%load('Processed_Batch_Data.mat');


% SWAPPED: Using ValData (Experiment 04) for training
Data = ValData;
t_bio = Data.Biomasse.t;
y_bio = Data.Biomasse.y;
t_glc = Data.Glucose.t;
y_glc = Data.Glucose.y;
var_bio = Data.Biomasse.var;
var_glc = Data.Glucose.var;

% 2. MODELL1 Modellkofiguration und Anfangswerte
% Modell1 expects concentrations (g/L), volume V is removed.
% Anfangswerte aus der jeweils ersten Messung.
x0 = [y_bio(1); y_glc(1)];

% Modell-Konfiguration
kinetic    = 3;      % 3 = Monod
%--------------------------------------------------------------------------
withOxygen = true;
%--------------------------------------------------------------------------
% Parametervektor [mu_max; K_S; Y_XS]
p0  = [0.3;  0.5;  0.15];    % Startwerte
pLB = [0.01; 0.01; 0.01];    % untere Schranken
pUB = [1.0;  500.0;  1.0];     % obere Schranken (K_S oben gelockert)

DOTstern = 0;
if withOxygen
    % Sauerstoffmessung (DOT)
    t_o2 = Data.O2.t;  
    y_o2 = Data.O2.y;  
    var_o2 = Data.O2.var;
    DOTstern = max(y_o2);

    % 3. Zustand ergaenzen: cO2 (Start-DOT aus O2-Messung zu Batch-Beginn)
    x0 = [x0; y_o2(1)];

    % Satrtwerte [Y_XO; KLa]
    %   Y_XO = Ertragskoeffizient Biomasse/O2
    %   KLa  = volumetrischer O2-Transferkoeffizient [1/h]
    p0  = [p0;  1.0; 200];
    pLB = [pLB; 0.01;   1];
    pUB = [pUB; 500.0; 1000];
end

% 3. MODELL1 Parameteridentifikation (WLS) -> uebung5.m / guete_pi_WLS.m

% Gütefunktion (Weighted Least Squares)
obj_fun = @(p) calculate_wls_error(p, x0, Data, kinetic, withOxygen);

% Durchfuehrung der Optimierung mit fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[p_opt, fval] = fmincon(obj_fun, p0, [], [], [], [], pLB, pUB, [], options);

% 4. MODELL1 Ausgabe und Speichern der Identifikationsergebnisse
fprintf('\n--- Parameteridentifikation abgeschlossen ---\n');
fprintf('Finaler WLS-Fehler: %.4f\n', fval);
fprintf('mu_max = %.4f 1/h\n', p_opt(1));
fprintf('K_S    = %.4f g/L\n', p_opt(2));
fprintf('Y_XS   = %.4f g/g\n', p_opt(3));
if withOxygen
    fprintf('Y_XO   = %.4f g/g\n', p_opt(4));
    fprintf('KLa    = %.4f 1/h\n', p_opt(5));
end

% Speichern der optimierten Parameter -> für main05
scriptDir = pwd;
saveDir = fullfile(scriptDir, '..', 'Daten');
savePath = fullfile(saveDir, 'p_opt.mat');
save(savePath,"p_opt");

%% 5. MODELL1 Visualisierung
% Simulation nur so lange laufen lassen, wie Batch-Daten vorhanden sind:
t_end       = max([t_bio(:)+1; t_glc(:)+1]);
t_sim       = linspace(0, t_end, 200);
options_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

[~, X_sim] = ode45(@(t, x) Modell1(t, x, p_opt, kinetic, withOxygen, DOTstern), t_sim, x0, options_ode);

c_X_sim   = X_sim(:, 1);
c_Glc_sim = X_sim(:, 2);
if withOxygen
    cO2_sim = X_sim(:, 3);
end

figure('Name', 'Task 4: Parameter Identification (Modell1)', 'Position', [200, 200, 900, 600]);
if withOxygen
    nRows = 3;
else
    nRows = 2;
end

subplot(nRows, 1, 1);
% Messung mit Fehlerbalken -> +/- 1 Standardabweichung aus der Messvarianz
errorbar(t_bio, y_bio, sqrt(var_bio), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim, c_X_sim, 'LineWidth', 2);
title('Biomass');
ylabel('c_X (g/L)');
legend("Messung \pm \sigma", "Simulation", Location="southeast");
xlim([0 9.5])
grid on;

subplot(nRows, 1, 2);
errorbar(t_glc, y_glc, sqrt(var_glc), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim, c_Glc_sim, 'LineWidth', 2);
title('Glucose');

ylabel('c_{Glc} (g/L)');
legend("Messung \pm \sigma", "Simulation", Location="northeast");
xlim([0 9.5])
grid on;

if withOxygen
    subplot(nRows, 1, 3);
    errorbar(t_o2, y_o2, sqrt(var_o2), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
    plot(t_sim, cO2_sim, 'LineWidth', 2);
    title('Sauerstoff (DOT)');
    ylabel('cO_2 (%)');
    legend("Messung \pm \sigma", "Simulation");
    xlim([0 9.5])
    grid on;
end
xlabel('BatchAge (h)');


%% --------------------------------------------------------------------------
%% MODELL 2 -> Visualisierung selber aufrufen
%% --------------------------------------------------------------------------
%projectRoot = pwd;
%load(fullfile(projectRoot,'..','Daten/Daten_Processed/Processed_Batch_Data.mat'));
%addpath(fullfile(projectRoot,'..','Modelle'),'-begin');

% SWAPPED: Using ValData (Experiment 04) for training
Data = ValData;

t_bio = Data.Biomasse.t;
y_bio = Data.Biomasse.y;
t_glc = Data.Glucose.t;
y_glc = Data.Glucose.y;
var_bio = Data.Biomasse.var;
var_glc = Data.Glucose.var;

t_am = Data.Ammonium.t;  y_am = Data.Ammonium.y;  var_am = Data.Ammonium.var;
t_ba = Data.Base.t;      y_ba = Data.Base.y;      var_ba = Data.Base.var;
t_o2 = Data.O2.t;        y_o2 = Data.O2.y;        var_o2 = Data.O2.var;
DOTstern = max(y_o2);


% 1. Anfangswerte und Parameter
% Zustaende: [cX; cGlc; cAm; cBase; cO2]
x0_m2 = [y_bio(1); y_glc(1); y_am(1); y_ba(1); y_o2(1)];

% Parameter: [mu_max; K_S; Y_XS; Y_Bam; Y_AmX; Y_XO; KLa]
p0_m2  = [0.3;  0.5;  0.15; 1.0;  0.05;  1.0; 200];
pLB_m2 = [0.01; 0.01; 0.01; 0.01; 0.001; 0.01;   1];
pUB_m2 = [1.0; 5.0; 1.0; 10; 1.0; 5.0; 1000]; 
%pUB KS von 5 auf 500 erweitert, weil er bei 5 und 50 an die obere Grenze gestoßen ist
%mu_max von 1 auf 10 erhöht


% Modell-Konfiguration
kinetic    = 3;      % 3 = Monod

% 2. Parameteridentifikation (WLS)
obj_fun_m2 = @(p) calculate_wls_error_m2(p, x0_m2, Data, kinetic);

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[p_opt_m2, fval_m2] = fmincon(obj_fun_m2, p0_m2, [], [], [], [], pLB_m2, pUB_m2, [], options);

% 3. Ausgabe und Speichern
fprintf('\n--- Modell2: Parameteridentifikation abgeschlossen ---\n');
fprintf('Finaler WLS-Fehler: %.4f\n', fval_m2);
fprintf('mu_max = %.4f 1/h\n',  p_opt_m2(1));
fprintf('K_S    = %.4f g/L\n',  p_opt_m2(2));
fprintf('Y_XS   = %.4f g/g\n',  p_opt_m2(3));
fprintf('Y_Bam  = %.4f\n',      p_opt_m2(4));
fprintf('Y_AmX  = %.4f\n',      p_opt_m2(5));
fprintf('Y_XO   = %.4f g/g\n',  p_opt_m2(6));
fprintf('KLa    = %.4f 1/h\n',  p_opt_m2(7));

% Speichern von Modell2
scriptDir = pwd;
saveDir = fullfile(scriptDir, '..', 'Daten');
save(fullfile(saveDir, 'p_opt_Modell2.mat'), 'p_opt_m2');

%% 4. Visualisierung Modell 2
t_end_m2 = max([t_bio(:)+1; t_glc(:)+1; t_am(:)+1; t_ba(:)+1; t_o2(:)+1]);
t_sim_m2 = linspace(0, t_end_m2, 300);

options_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[~, X2] = ode15s(@(t, x) Modell2(t, x, p_opt_m2, kinetic, DOTstern), t_sim_m2, x0_m2, options_ode);

figure('Name', 'Task 4: Parameter Identification (Modell2: +Base +O2)', 'Position', [250, 80, 950, 950]);

subplot(5, 1, 1);
errorbar(t_bio, y_bio, sqrt(var_bio), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:, 1), 'LineWidth', 2);
title('Modell2 - Biomasse');
ylabel('c_X (g/L)');
legend("Messung \pm \sigma", "Simulation");
xlim([0 10])
grid on;

subplot(5, 1, 2);
errorbar(t_glc, y_glc, sqrt(var_glc), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:, 2), 'LineWidth', 2);
title('Modell2 - Glucose');
ylabel('c_{Glc} (g/L)');
legend("Messung \pm \sigma", "Simulation");
xlim([0 10])
grid on;

subplot(5, 1, 3);
errorbar(t_am, y_am, sqrt(var_am), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:, 3), 'LineWidth', 2);
title('Modell2 - Ammonium');
ylabel('c_{Am} (g/L)');
legend("Messung \pm \sigma", "Simulation");
xlim([0 10])
grid on;

subplot(5, 1, 4);
errorbar(t_ba, y_ba, sqrt(var_ba), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:, 4), 'LineWidth', 2);
title('Modell2 - Base');
ylabel('c_{Base}');
legend("Messung \pm \sigma", "Simulation");
xlim([0 10])
grid on;

subplot(5, 1, 5);
errorbar(t_o2, y_o2, sqrt(var_o2), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 4); hold on;
plot(t_sim_m2, X2(:, 5), 'LineWidth', 2);
title('Modell2 - Sauerstoff (DOT)');
xlabel('BatchAge (h)');
ylabel('cO_2 (%)');
legend("Messung \pm \sigma", "Simulation");
xlim([0 10])
grid on;




%% Local Objective Function -> guete_pi_WLS (WLS)

% function J = calculate_wls_error(p, x0, t_mes, y_bio, y_glc, var_bio, var_glc, kinetic, withOxygen)
%     % Erzwingt Auswertung von ode15s exakt an den Messzeitpunkten.
%     try
%         [~, X_sim] = ode15s(@(t, x) Modell1(t, x, p, kinetic, withOxygen), t_mes, x0);
%         c_X_sim   = X_sim(:, 1);
%         c_Glc_sim = X_sim(:, 2);
% 
%         % Sicherstellen, dass alle Groessen Spaltenvektoren sind
%         y_bio = y_bio(:); y_glc = y_glc(:);
%         var_bio = var_bio(:); var_glc = var_glc(:);
% 
%         err_bio = sum(((y_bio - c_X_sim).^2) ./ var_bio);
%         err_glc = sum(((y_glc - c_Glc_sim).^2) ./ var_glc);
% 
%         J = err_bio + err_glc;
%     catch
%         J = 1e6;   % Strafterm bei Integrationsfehler
%     end
% end
%% Local Objective Function -> guete_pi_WLS (WLS) -- GEGENCHECKEN!
% Hinweis: Die O2-Messung (DOT) hat einen viel feineren Zeitvektor als die
% Offline-Groessen Biomasse/Glucose. Die Guetefunktion simuliert daher
% einmalig ueber das vereinigte Zeitraster aller Messungen und ordnet die
% Simulationswerte anschliessend den jeweiligen Messzeitpunkten zu.

function J = calculate_wls_error(p, x0, Data, kinetic, withOxygen)
    % Messgroessen (jede Groesse mit eigenem Zeitvektor)
    t_bio = Data.Biomasse.t(:); y_bio = Data.Biomasse.y(:); var_bio = Data.Biomasse.var(:);
    t_glc = Data.Glucose.t(:);  y_glc = Data.Glucose.y(:);  var_glc = Data.Glucose.var(:);
    DOTstern = 0;
    if withOxygen
        t_o2 = Data.O2.t(:); y_o2 = Data.O2.y(:); var_o2 = Data.O2.var(:);
        t_all = unique([t_bio; t_glc; t_o2]);   % vereinigtes Zeitraster
        DOTstern = max(y_o2);
    else
        t_all = unique([t_bio; t_glc]);
    end

    % Einmalige Simulation ueber das vereinigte Zeitraster.
    % (ode15s liefert bei >= 3 Zeitpunkten die Loesung exakt an diesen Punkten.)
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    try
        [~, X_sim] = ode15s(@(t, x) Modell1(t, x, p, kinetic, withOxygen, DOTstern), t_all, x0, options);
    catch ME
        warning('Modell1-Integration fehlgeschlagen: %s', ME.message);
        J = 1e6;
        return;
    end


    % Simulationswerte den jeweiligen Messzeitpunkten zuordnen
    [~, iBio] = ismember(t_bio, t_all);
    [~, iGlc] = ismember(t_glc, t_all);
    % 
    % err_bio = sum(((y_bio - X_sim(iBio, 1)).^2) ./ var_bio);
    % err_glc = sum(((y_glc - X_sim(iGlc, 2)).^2) ./ var_glc);
    % J = err_bio + err_glc;
    % 
    % if withOxygen
    %     [~, iO2] = ismember(t_o2, t_all);
    %     err_o2 = sum(((y_o2 - X_sim(iO2, 3)).^2) ./ var_o2);
    %     J = J +  err_o2;
    % end

    err_bio = sum(((y_bio - X_sim(iBio, 1)).^2) ./ var_bio.^2);
    err_glc = sum(((y_glc - X_sim(iGlc, 2)).^2) ./ var_glc.^2);
    J = err_bio + err_glc;

    if withOxygen
        [~, iO2] = ismember(t_o2, t_all);
        err_o2 = 1/20 .* sum(((y_o2 - X_sim(iO2, 3)).^2) ./ var_o2.^2);
        J = J +  err_o2;
    end

    if ~isfinite(J); J = 1e6; end
end


function J = calculate_wls_error_m2(p, x0, Data, kinetic)
    % Modell2: Biomasse, Glucose, Ammonium, Base, O2.
    t_bio = Data.Biomasse.t(:); y_bio = Data.Biomasse.y(:); var_bio = Data.Biomasse.var(:);
    t_glc = Data.Glucose.t(:);  y_glc = Data.Glucose.y(:);  var_glc = Data.Glucose.var(:);
    t_am  = Data.Ammonium.t(:); y_am  = Data.Ammonium.y(:); var_am  = Data.Ammonium.var(:);
    t_ba  = Data.Base.t(:);     y_ba  = Data.Base.y(:);     var_ba  = Data.Base.var(:);
    t_o2  = Data.O2.t(:);       y_o2  = Data.O2.y(:);       var_o2  = Data.O2.var(:);

    t_all = unique([t_bio; t_glc; t_am; t_ba; t_o2]);   % vereinigtes Zeitraster
    DOTstern = max(y_o2);

    % Einmalige Simulation ueber das vereinigte Zeitraster.
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    try
        [~, X_sim] = ode15s(@(t, x) Modell2(t, x, p, kinetic, DOTstern), t_all, x0, options);
    catch ME
        warning('Modell2-Integration fehlgeschlagen: %s', ME.message);
        J = 1e6;
        return;
    end

    % Simulationswerte den jeweiligen Messzeitpunkten zuordnen
    [~, iBio] = ismember(t_bio, t_all);
    [~, iGlc] = ismember(t_glc, t_all);
    [~, iAm]  = ismember(t_am,  t_all);
    [~, iBa]  = ismember(t_ba,  t_all);
    [~, iO2]  = ismember(t_o2,  t_all);

    % err_bio = sum(((y_bio - X_sim(iBio, 1)).^2) ./ var_bio);
    % err_glc = sum(((y_glc - X_sim(iGlc, 2)).^2) ./ var_glc);
    % err_am  = sum(((y_am  - X_sim(iAm,  3)).^2) ./ var_am);
    % err_ba  = sum(((y_ba  - X_sim(iBa,  4)).^2) ./ var_ba);
    % err_o2  = sum(((y_o2  - X_sim(iO2,  5)).^2) ./ var_o2);
    % J = err_bio + err_glc + err_am + err_ba + err_o2;

    err_bio = sum(((y_bio - X_sim(iBio, 1)).^2) ./ var_bio.^2);
    err_glc = sum(((y_glc - X_sim(iGlc, 2)).^2) ./ var_glc.^2);
    err_am  = sum(((y_am  - X_sim(iAm,  3)).^2) ./ var_am.^2);
    err_ba  = sum(((y_ba  - X_sim(iBa,  4)).^2) ./ var_ba.^2);
    err_o2  = 1/20 .* sum(((y_o2  - X_sim(iO2,  5)).^2) ./ var_o2.^2);
    J = err_bio + err_glc + err_am + err_ba + err_o2;

    if ~isfinite(J); J = 1e6; end
end
