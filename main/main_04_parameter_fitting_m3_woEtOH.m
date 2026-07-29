%% main_04_parameter_fitting_m3_woEtOH.m
%  Parameterfitting Modell 3 (Fed-Batch, ohne Ethanol) auf MEHREREN
%  Trainingsexperimenten (gemeinsamer Parametersatz) + Validierung auf
%  einem gehaltenen Experiment. Basisvariante (ein fmincon-Start).
%
%  Erwartet aus main_01_data_preprocessing_m3.m:
%    TrainData : 1xN Cell-Array (je ein Experiment-Struct)
%    ValData   : ein Experiment-Struct
% =========================================================================
clear; clc; close all;

projectRoot = pwd;
load(fullfile(projectRoot, '..', 'Daten', 'Daten_Processed', 'Processed_FedBatch_Modell3.mat'));
addpath(fullfile(projectRoot, '..', 'Modelle'), '-begin');
rehash; clear Modell3_woEtOH

% Parameter [mumax KS YXS YAmX YPhX YB_Am KLa YXO]
namen = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO'};
p0  = [0.30, 0.50, 0.15, 0.05,  0.02,  1.0, 200, 1.0];
pLB = [0.01, 0.01, 0.05, 0.001, 0.001, 0.1,  10, 0.1];
pUB = [1.00, 5.00, 1.00, 1.000, 1.000, 5.0, 800, 5.0];

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', 5000, ...
                       'FiniteDifferenceType', 'central', ...
                       'FiniteDifferenceStepSize', 1e-6, ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance', 1e-10);

obj_fun = @(p) wls_m3_multi(p, TrainData);
[p_opt, fval] = fmincon(obj_fun, p0, [], [], [], [], pLB, pUB, [], options);

fprintf('\n--- Modell3_woEtOH: Parameteridentifikation abgeschlossen ---\n');
fprintf('Training %s | Validierung %d\n', mat2str(trainNums), valNum);
fprintf('Finaler WLS-Fehler (Training):    %.4f\n', fval);
for i = 1:numel(p_opt)
    fprintf('%-8s = %.4f\n', namen{i}, p_opt(i));
end
fprintf('WLS-Fehler (Validierung %d):      %.4f\n', valNum, wls_m3_single(p_opt, ValData));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt_Modell3_woEtOH.mat'), 'p_opt');

% Visualisierung: jedes Trainingsexperiment + Validierung
for e = 1:numel(TrainData)
    plot_fit_m3(TrainData{e}, p_opt, sprintf('Modell3 | Training %s', TrainData{e}.Name));
end
plot_fit_m3(ValData, p_opt, sprintf('Modell3 | Validierung %s', ValData.Name));


%% ======================================================================
%% WLS-Guetefunktionen
%% ======================================================================
function J = wls_m3_multi(p, Train)
% Summe der WLS-Fehler ueber alle Trainingsexperimente (gemeinsames p).
    J = 0;
    for e = 1:numel(Train)
        J = J + wls_m3_single(p, Train{e});
    end
    if ~isfinite(J), J = 1e8; end
end

function J = wls_m3_single(p, D)
% WLS eines Experiments fuer Modell3_woEtOH.
% Zustaende: [V; mX; mGlc; mAm; mPh; mB; DOT].
% Integration ab x0 bei t0 = u(1,1); bewertet werden die in D vorhandenen
% Messpunkte (bei gefensterten Daten also nur die behaltenen).
    u = D.u;  x0 = D.x0;  DOTstern = D.DOTstern;

    M = { D.Biomasse, 'Biomasse', 2, true;  ...
          D.Glucose,  'Glucose',  3, true;  ...
          D.Ammonium, 'Ammonium', 4, true;  ...
          D.Phosphat, 'Phosphat', 5, true;  ...
          D.Base,     'Base',     6, false; ...
          D.O2,       'DOT',      7, false  };
    wmode = 'mean';
    wsig  = [1, 1, 1, 1, 1, 1];

    % Vereinigtes Zeitraster inkl. Startzeit t0 und Feed-Sprungzeiten
    t0    = u(1,1);
    t_all = t0;
    for i = 1:size(M,1), t_all = [t_all; M{i,1}.t(:)]; end
    tu = u(1,:).';
    tu = tu(tu >= min(t_all) & tu <= max(t_all));
    t_all = unique([t_all; tu]);
    if t_all(1) > t0, t_all = [t0; t_all]; end

    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);
    try
        [~, X] = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p,DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:))), J = 1e8; return; end

    V = X(:,1);
    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1}; name = M{i,2}; idxS = M{i,3}; divByV = M{i,4};
        if isempty(mess.t), continue; end
        [tf, iT] = ismember(mess.t(:), t_all);
        if any(~tf), warning('Messzeitpunkt fuer %s nicht gefunden.', name); J = 1e8; return; end
        y_sim = X(iT, idxS);
        if divByV, y_sim = y_sim ./ V(iT); end
        r = (mess.y(:) - y_sim) ./ sqrt(max(mess.var(:), eps));
        switch wmode
            case 'sum',  contrib = sum(r.^2);
            case 'mean', contrib = mean(r.^2);
        end
        J = J + wsig(i) * contrib;
    end
    if ~isfinite(J), J = 1e8; end
end


%% ======================================================================
%% Visualisierung
%% ======================================================================
function plot_fit_m3(D, p, titel)
    u = D.u;  x0 = D.x0;  DOTstern = D.DOTstern;

    t_start = u(1,1);
    t_end   = max(D.Biomasse.t(:)) + 1;
    t_sim   = linspace(t_start, t_end, 300);
    [~, X]  = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p,DOTstern), ...
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
