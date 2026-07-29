%% main_02_simulation_m3_woEtOH.m
%  Simulations-Check Modell 3 (Fed-Batch, ohne Ethanol) VOR dem Fitting.
%  Simuliert mit den Startschaetzwerten (p_guess) und stellt fuer JEDES
%  Experiment (alle Trainingsdatensaetze + Validierung) Messung vs.
%  Simulation gegenueber.
%
%  Erwartet aus main_01_data_preprocessing_m3.m:
%    TrainData : 1xN Cell-Array (je ein Experiment-Struct)
%    ValData   : ein Experiment-Struct
% =========================================================================
clear; clc; close all;

scriptDir = pwd;
loadPath  = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed', 'Processed_FedBatch_Modell3.mat');
addpath(fullfile(scriptDir, '..', 'Modelle'), '-begin');
rehash; clear Modell3_woEtOH

if ~isfile(loadPath)
    error('Processed_FedBatch_Modell3.mat nicht gefunden. Zuerst main_01_..._m3 ausfuehren.');
end
load(loadPath);
fprintf('Daten geladen:\n  %s\n', loadPath);

% Startschaetzwerte [mumax,KS,YXS,YAmX,YPhX,YB_Am,KLa,YXO]
p_guess = [0.3, 0.5, 0.15, 0.05, 0.02, 1.0, 200, 1.0];

% Alle Datensaetze (Training + Validierung)
alle  = [TrainData, {ValData}];
rolle = [repmat({'Training'}, 1, numel(TrainData)), {'Validierung'}];

for e = 1:numel(alle)
    D   = alle{e};
    tag = sprintf('%s %s', rolle{e}, D.Name);
    sim_plot_m3(D, p_guess, sprintf('Modell3 | %s', tag));
end

fprintf('[OK] Simulations-Check fuer %d Datensaetze erstellt.\n', numel(alle));


%% ======================================================================
%% Simulation + Plot
%% ======================================================================
function sim_plot_m3(D, p, titel)
    u  = D.u;  x0 = D.x0;  DOTstern = D.DOTstern;

    t_start = u(1,1);
    t_end   = max(D.Biomasse.t(:)) + 1;
    t_sim   = linspace(t_start, t_end, 300);

    [~, X] = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p,DOTstern), ...
                    t_sim, x0, odeset('RelTol',1e-5,'AbsTol',1e-7));

    V   = X(:,1);
    cX  = X(:,2)./V;  cGlc = X(:,3)./V;  cAm = X(:,4)./V;
    cPh = X(:,5)./V;  mB   = X(:,6);     DOT = X(:,7);

    figure('Name', titel, 'Position', [200 40 950 1000]);
    plot_row(1, D.Biomasse, t_sim, cX,   [titel ' - Biomasse'], 'c_X (g/L)');
    plot_row(2, D.Glucose,  t_sim, cGlc, 'Glucose',   'c_{Glc} (g/L)');
    plot_row(3, D.Ammonium, t_sim, cAm,  'Ammonium',  'c_{Am} (g/L)');
    plot_row(4, D.Phosphat, t_sim, cPh,  'Phosphat',  'c_{Ph} (g/L)');
    plot_row(5, D.Base,     t_sim, mB,   'Base',      'm_B');
    plot_row(6, D.O2,       t_sim, DOT,  'DOT',       'DOT (%)');
    xlabel('BatchAge (h)');
end


function plot_row(row, mess, t_sim, y_sim, name, ylab)
    subplot(6,1,row);
    errorbar(mess.t, mess.y, sqrt(mess.var), 'o', 'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, y_sim, 'LineWidth', 2);
    title(name); ylabel(ylab);
    legend('Messung \pm \sigma','Simulation','Location','best'); grid on;
end
