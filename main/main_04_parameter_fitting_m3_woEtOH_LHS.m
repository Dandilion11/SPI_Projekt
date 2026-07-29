%% main_04_parameter_fitting_m3_woEtOH_LHS.m
%  Parameterfitting Modell 3 (Fed-Batch, ohne Ethanol) mit LHS-Multistart,
%  jetzt mit MULTI-DATENSATZ-TRAINING (gemeinsamer Parametersatz) +
%  Validierung auf einem gehaltenen Experiment.
%
%  Beibehaltenes Feature der bisherigen LHS-Variante ("from feed"):
%   - Pro Experiment wird der erste Feed-Zeitpunkt bestimmt und nur ab
%     t >= t_feed (bis t_hi) gefittet. Die ODE integriert weiterhin ab x0
%     bei t0 = u(1,1) -> Anfangsbedingung bleibt korrekt, nur die dichten,
%     flachen Baseline-Punkte vor dem Feed werden nicht bewertet.
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

% ---- Fenster- und LHS-Konfiguration ----
t_hi  = inf;   % obere Fenstergrenze (inf = kein Schwanz-Cut; z.B. 60 zum Kappen)
N_lhs = 100;   % LHS-Screening-Punkte (zum Testen z.B. 30)
K_opt = 3;     % beste Startpunkte -> teurer fmincon (zum Testen z.B. 2)

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

% Pro Trainingsexperiment ab dem eigenen Feed-Start fenstern
TrainFit = cell(size(TrainData));
for e = 1:numel(TrainData)
    t_feed      = feed_start(TrainData{e}.u);
    TrainFit{e} = window_data(TrainData{e}, t_feed, t_hi);
    fprintf('%s: Feed-Start t_feed = %.3f h -> Fit ab hier.\n', TrainData{e}.Name, t_feed);
    report_window(TrainData{e}, TrainFit{e});
end

% Validierung ebenfalls fenstern (nur fuer den WLS-Vergleich)
t_feed_val = feed_start(ValData.u);
ValFit     = window_data(ValData, t_feed_val, t_hi);

obj_fun = @(p) wls_m3_multi(p, TrainFit);
[p_opt, fval] = lhs_multistart(obj_fun, p0, pLB, pUB, options, N_lhs, K_opt);

fprintf('\n--- Modell3_woEtOH (LHS, from feed): Fitting abgeschlossen ---\n');
fprintf('Training %s | Validierung %d\n', mat2str(trainNums), valNum);
fprintf('Finaler WLS-Fehler (Training, ab Feed): %.4f\n', fval);
for i = 1:numel(p_opt)
    fprintf('%-8s = %.4f\n', namen{i}, p_opt(i));
end
fprintf('WLS-Fehler (Validierung %d, ab Feed):   %.4f\n', valNum, wls_m3_single(p_opt, ValFit));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt_Modell3_woEtOH.mat'), 'p_opt');

% Visualisierung ueber den VOLLEN Horizont, mit Feed-Linie
for e = 1:numel(TrainData)
    plot_fit_m3(TrainData{e}, p_opt, feed_start(TrainData{e}.u), t_hi, ...
                sprintf('Modell3 LHS | Training %s', TrainData{e}.Name));
end
plot_fit_m3(ValData, p_opt, t_feed_val, t_hi, ...
            sprintf('Modell3 LHS | Validierung %s', ValData.Name));


%% ======================================================================
%% LHS-Multistart
%% ======================================================================
function [p_opt, fval] = lhs_multistart(obj_fun, p0, pLB, pUB, options, N_lhs, K_opt)
    wasCol = iscolumn(p0);
    p0r  = p0(:).';  pLBr = pLB(:).';  pUBr = pUB(:).';
    d    = numel(p0r);

    L     = lhsdesign(N_lhs, d);
    logLB = log10(pLBr);  logUB = log10(pUBr);
    P     = 10.^(logLB + L .* (logUB - logLB));

    fprintf('LHS-Screening ueber %d Punkte ...\n', N_lhs);
    Jscreen = inf(N_lhs,1);
    for k = 1:N_lhs
        try, Jscreen(k) = obj_fun(P(k,:)); catch, end
    end
    [~, order] = sort(Jscreen);
    K      = min(K_opt, N_lhs);
    Pstart = [p0r; P(order(1:K), :)];

    p_opt = p0r;  fval = inf;
    for k = 1:size(Pstart,1)
        fprintf('--- fmincon Start %d/%d ---\n', k, size(Pstart,1));
        try
            [pk, Jk] = fmincon(obj_fun, Pstart(k,:), [], [], [], [], pLBr, pUBr, [], options);
            if Jk < fval, fval = Jk;  p_opt = pk; end
        catch
        end
    end
    if wasCol, p_opt = p_opt(:); end
end


%% ======================================================================
%% Feed-Start + Fensterung
%% ======================================================================
function t_feed = feed_start(u)
% Erster Zeitpunkt, zu dem eine Substrat-Feedrate > 0 wird.
% Feedraten in u: Zeile 2 = uAm, 4 = uPh, 6 = uGlc.
    feedRows = [2 4 6];
    feedOn   = any(u(feedRows,:) > 0, 1);
    idx      = find(feedOn, 1, 'first');
    if isempty(idx), t_feed = u(1,1); else, t_feed = u(1, idx); end
end

function D = window_data(D, t_lo, t_hi)
% Behaelt aus jeder Messgroesse nur die Punkte mit t_lo <= t <= t_hi.
% u, x0, DOTstern, Name, V bleiben unveraendert.
    flds = {'Biomasse','Glucose','Ammonium','Phosphat','Base','O2'};
    for i = 1:numel(flds)
        f = flds{i};  m = D.(f);
        keep = (m.t >= t_lo) & (m.t <= t_hi);
        m.t = m.t(keep);  m.y = m.y(keep);  m.var = m.var(keep);
        D.(f) = m;
    end
end

function report_window(Dfull, Dcut)
    flds = {'Biomasse','Glucose','Ammonium','Phosphat','Base','O2'};
    for i = 1:numel(flds)
        f = flds{i};
        fprintf('   %-9s: %d -> %d\n', f, numel(Dfull.(f).t), numel(Dcut.(f).t));
    end
end


%% ======================================================================
%% WLS-Guetefunktionen
%% ======================================================================
function J = wls_m3_multi(p, Train)
    J = 0;
    for e = 1:numel(Train)
        J = J + wls_m3_single(p, Train{e});
    end
    if ~isfinite(J), J = 1e8; end
end

function J = wls_m3_single(p, D)
% WLS eines (ggf. gefensterten) Experiments fuer Modell3_woEtOH.
% Integration ab x0 bei t0 = u(1,1); bewertet werden nur die in D
% vorhandenen (behaltenen) Messpunkte.
    u = D.u;  x0 = D.x0;  DOTstern = D.DOTstern;

    M = { D.Biomasse, 'Biomasse', 2, true;  ...
          D.Glucose,  'Glucose',  3, true;  ...
          D.Ammonium, 'Ammonium', 4, true;  ...
          D.Phosphat, 'Phosphat', 5, true;  ...
          D.Base,     'Base',     6, false; ...
          D.O2,       'DOT',      7, false  };
    wmode = 'mean';
    wsig  = [1, 1, 1, 1, 1, 1];

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
%% Visualisierung (voller Horizont, Feed-Linie)
%% ======================================================================
function plot_fit_m3(D, p, t_feed, t_hi, titel)
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
    plot_row(1, D.Biomasse, t_sim, cX,   [titel ' - Biomasse'], 'c_X (g/L)',     t_feed, t_hi);
    plot_row(2, D.Glucose,  t_sim, cGlc, 'Glucose',   'c_{Glc} (g/L)',           t_feed, t_hi);
    plot_row(3, D.Ammonium, t_sim, cAm,  'Ammonium',  'c_{Am} (g/L)',            t_feed, t_hi);
    plot_row(4, D.Phosphat, t_sim, cPh,  'Phosphat',  'c_{Ph} (g/L)',            t_feed, t_hi);
    plot_row(5, D.Base,     t_sim, mB,   'Base',      'm_B',                     t_feed, t_hi);
    plot_row(6, D.O2,       t_sim, DOT,  'DOT',       'DOT (%)',                 t_feed, t_hi);
    xlabel('BatchAge (h)');
end

function plot_row(row, mess, t_sim, y_sim, name, ylab, t_feed, t_hi)
    subplot(6,1,row);
    errorbar(mess.t, mess.y, sqrt(mess.var), 'o', 'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, y_sim, 'LineWidth', 2);
    xline(t_feed, '--k', 'Feed', 'LabelVerticalAlignment','middle');
    if nargin >= 8 && isfinite(t_hi)
        xline(t_hi, '--k', 't_{hi}', 'LabelVerticalAlignment','top');
    end
    title(name); ylabel(ylab);
    legend('Messung \pm \sigma','Simulation','Location','best'); grid on;
end
