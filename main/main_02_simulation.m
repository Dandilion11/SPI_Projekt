%% main_02_simulation.m
%  Simulations-Check fuer Modell 1 & 2 VOR dem Fitting.
%  Simuliert mit den Startschaetzwerten (p_guess) und stellt fuer JEDES
%  Experiment (alle Trainingsdatensaetze + Validierung) Messung vs.
%  Simulation gegenueber.
%
%  Erwartet aus main_01_data_preprocessing.m:
%    TrainData : 1xN Cell-Array (je ein Experiment-Struct)
%    ValData   : ein Experiment-Struct
% =========================================================================
clear; clc; close all;

scriptDir = pwd;
loadPath  = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat');
addpath(fullfile(scriptDir, '..', 'Modelle'), '-begin');
rehash; clear Modell1 Modell2

if ~isfile(loadPath)
    error('Processed_Batch_Data.mat nicht gefunden. Zuerst main_01 ausfuehren.');
end
load(loadPath);
fprintf('Daten geladen:\n  %s\n', loadPath);

kinetic    = 3;      % 3 = Monod
withOxygen = true;

% Startschaetzwerte (die echten Werte kommen aus main_04)
p_guess_m1 = [0.3, 0.5, 0.15];              % [mu_max, K_S, Y_XS]
if withOxygen
    p_guess_m1 = [p_guess_m1, 1.0, 200];    % + [Y_XO, KLa]
end
p_guess_m2 = [0.3, 0.5, 0.15, 1.0, 0.05, 1.0, 200]; % [mu_max,K_S,Y_XS,Y_Bam,Y_AmX,Y_XO,KLa]

% Alle Datensaetze zusammenstellen (Training + Validierung)
alle  = [TrainData, {ValData}];
rolle = [repmat({'Training'}, 1, numel(TrainData)), {'Validierung'}];

% Fuer jedes Experiment Modell 1 und Modell 2 simulieren und plotten
for e = 1:numel(alle)
    D   = alle{e};
    tag = sprintf('%s %s', rolle{e}, D.Name);
    sim_plot_m1(D, p_guess_m1, kinetic, withOxygen, sprintf('Modell1 | %s', tag));
    sim_plot_m2(D, p_guess_m2, kinetic,             sprintf('Modell2 | %s', tag));
end

fprintf('[OK] Simulations-Check fuer %d Datensaetze erstellt (Modell 1 & 2).\n', numel(alle));


%% ======================================================================
%% Simulation + Plot
%% ======================================================================
function sim_plot_m1(D, p, kinetic, withOxygen, titel)
    x0 = [D.Biomasse.y(1); D.Glucose.y(1)];
    DOTstern = 0;
    if withOxygen, x0 = [x0; D.O2.y(1)]; DOTstern = max(D.O2.y); end

    t_end = max([D.Biomasse.t(:); D.Glucose.t(:)]) + 1;
    t_sim = linspace(0, t_end, 200);
    [~, X] = ode15s(@(t,x) Modell1(t,x,p,kinetic,withOxygen,DOTstern), ...
                    t_sim, x0, odeset('RelTol',1e-5,'AbsTol',1e-7));

    n = 2 + double(withOxygen);
    figure('Name', titel, 'Position', [200 150 850 600]);
    subplot(n,1,1);
    errorbar(D.Biomasse.t, D.Biomasse.y, sqrt(D.Biomasse.var), 'o', 'MarkerSize', 4); hold on;
    plot(t_sim, X(:,1), 'LineWidth', 2); title([titel ' - Biomasse']);
    ylabel('c_X (g/L)'); legend('Messung \pm \sigma','Simulation'); grid on;
    subplot(n,1,2);
    errorbar(D.Glucose.t, D.Glucose.y, sqrt(D.Glucose.var), 'o', 'MarkerSize', 4); hold on;
    plot(t_sim, X(:,2), 'LineWidth', 2); title('Glucose');
    ylabel('c_{Glc} (g/L)'); legend('Messung \pm \sigma','Simulation'); grid on;
    if withOxygen
        subplot(n,1,3);
        errorbar(D.O2.t, D.O2.y, sqrt(D.O2.var), 'o', 'MarkerSize', 4); hold on;
        plot(t_sim, X(:,3), 'LineWidth', 2); title('Sauerstoff (DOT)');
        ylabel('cO_2 (%)'); legend('Messung \pm \sigma','Simulation'); grid on;
    end
    xlabel('BatchAge (h)');
end


function sim_plot_m2(D, p, kinetic, titel)
    x0 = [D.Biomasse.y(1); D.Glucose.y(1); D.Ammonium.y(1); D.Base.y(1); D.O2.y(1)];
    DOTstern = max(D.O2.y);

    t_end = max([D.Biomasse.t(:); D.Glucose.t(:); D.Ammonium.t(:); D.Base.t(:); D.O2.t(:)]) + 1;
    t_sim = linspace(0, t_end, 300);
    [~, X] = ode15s(@(t,x) Modell2(t,x,p,kinetic,DOTstern), ...
                    t_sim, x0, odeset('RelTol',1e-5,'AbsTol',1e-7));

    kanal = {D.Biomasse,1,'Biomasse','c_X (g/L)'; ...
             D.Glucose, 2,'Glucose', 'c_{Glc} (g/L)'; ...
             D.Ammonium,3,'Ammonium','c_{Am} (g/L)'; ...
             D.Base,    4,'Base',    'c_{Base}'; ...
             D.O2,      5,'DOT',     'cO_2 (%)'};
    figure('Name', titel, 'Position', [250 60 900 900]);
    for k = 1:5
        subplot(5,1,k);
        errorbar(kanal{k,1}.t, kanal{k,1}.y, sqrt(kanal{k,1}.var), 'o', 'MarkerSize', 4); hold on;
        plot(t_sim, X(:, kanal{k,2}), 'LineWidth', 2);
        title([titel ' - ' kanal{k,3}]); ylabel(kanal{k,4});
        legend('Messung \pm \sigma','Simulation'); grid on;
    end
    xlabel('BatchAge (h)');
end
