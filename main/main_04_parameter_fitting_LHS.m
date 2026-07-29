%% main_04_parameter_fitting_LHS.m
%  Parameterfitting Modell 1 & 2 mit LHS-Multistart, jetzt mit
%  MULTI-DATENSATZ-TRAINING (gemeinsamer Parametersatz ueber mehrere
%  Experimente) + Validierung auf einem gehaltenen Experiment.
%
%  Beibehaltene Features aus der bisherigen LHS-Variante:
%   - LHS-Multistart um fmincon (billiges Screening -> beste Startpunkte)
%   - nCutDOT: fuehrende DOT-Punkte werden NICHT gefittet (x0 bleibt korrekt,
%     Startwert DOT = erster behaltener Punkt). Pro Experiment angewandt.
%
%  Erwartet aus main_01_data_preprocessing.m:
%    TrainData : 1xN Cell-Array (je ein Experiment-Struct)
%    ValData   : ein Experiment-Struct
% =========================================================================
clear; clc; close all;

projectRoot = pwd;
load(fullfile(projectRoot, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat'));
addpath(fullfile(projectRoot, '..', 'Modelle'), '-begin');
rehash; clear Modell1 Modell2

kinetic = 3;   % 3 = Monod

% ---- LHS-Multistart-Konfiguration ----
N_lhs = 100;   % LHS-Screening-Punkte (billig; zum Testen z.B. 30)
K_opt = 3;     % beste Startpunkte -> teurer fmincon (zum Testen z.B. 2)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

%% ======================================================================
%% MODELL 1
%% ======================================================================
withOxygen = true;
nCutDOT    = 1;    % Anzahl fuehrender DOT-Punkte, die NICHT gefittet werden

% Parametervektor [mu_max; K_S; Y_XS] (+ [Y_XO; KLa] bei withOxygen)
p0  = [0.3;  0.5;  0.15];
pLB = [0.01; 0.01; 0.01];
pUB = [1.0;  500;  1.0];
if withOxygen
    p0  = [p0;  1.0; 200];
    pLB = [pLB; 0.01;  1];
    pUB = [pUB; 500; 1000];
end

obj_m1 = @(p) wls_m1_multi(p, TrainData, kinetic, withOxygen, nCutDOT);
[p_opt, fval] = lhs_multistart(obj_m1, p0, pLB, pUB, options, N_lhs, K_opt);

fprintf('\n--- Modell1 (LHS): Parameteridentifikation abgeschlossen ---\n');
fprintf('Training %s | Validierung %d\n', mat2str(trainNums), valNum);
fprintf('Finaler WLS-Fehler (Training):    %.4f\n', fval);
fprintf('mu_max = %.4f 1/h\n', p_opt(1));
fprintf('K_S    = %.4f g/L\n', p_opt(2));
fprintf('Y_XS   = %.4f g/g\n', p_opt(3));
if withOxygen
    fprintf('Y_XO   = %.4f g/g\n', p_opt(4));
    fprintf('KLa    = %.4f 1/h\n', p_opt(5));
end
fprintf('WLS-Fehler (Validierung %d):      %.4f\n', valNum, ...
        wls_m1_single(p_opt, ValData, kinetic, withOxygen, nCutDOT));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt.mat'), 'p_opt');

for e = 1:numel(TrainData)
    plot_fit_m1(TrainData{e}, p_opt, kinetic, withOxygen, nCutDOT, ...
                sprintf('Modell1 LHS | Training %s', TrainData{e}.Name));
end
plot_fit_m1(ValData, p_opt, kinetic, withOxygen, nCutDOT, ...
            sprintf('Modell1 LHS | Validierung %s', ValData.Name));


%% ======================================================================
%% MODELL 2
%% ======================================================================
nCutDOT_m2 = 1;   % fuehrende DOT-Punkte, die NICHT gefittet werden

% Parameter [mu_max; K_S; Y_XS; Y_Bam; Y_AmX; Y_XO; KLa]
p0_m2  = [0.3;  0.5;  0.15; 1.0;  0.05;  1.0; 200];
pLB_m2 = [0.01; 0.01; 0.01; 0.01; 0.001; 0.01;   1];
pUB_m2 = [1.0;  5.0;  1.0;  10;   1.0;   100; 1000];

obj_m2 = @(p) wls_m2_multi(p, TrainData, kinetic, nCutDOT_m2);
[p_opt_m2, fval_m2] = lhs_multistart(obj_m2, p0_m2, pLB_m2, pUB_m2, options, N_lhs, K_opt);

fprintf('\n--- Modell2 (LHS): Parameteridentifikation abgeschlossen ---\n');
fprintf('Training %s | Validierung %d\n', mat2str(trainNums), valNum);
fprintf('Finaler WLS-Fehler (Training):    %.4f\n', fval_m2);
fprintf('mu_max = %.4f 1/h\n', p_opt_m2(1));
fprintf('K_S    = %.4f g/L\n', p_opt_m2(2));
fprintf('Y_XS   = %.4f g/g\n', p_opt_m2(3));
fprintf('Y_Bam  = %.4f\n',     p_opt_m2(4));
fprintf('Y_AmX  = %.4f\n',     p_opt_m2(5));
fprintf('Y_XO   = %.4f g/g\n', p_opt_m2(6));
fprintf('KLa    = %.4f 1/h\n', p_opt_m2(7));
fprintf('WLS-Fehler (Validierung %d):      %.4f\n', valNum, ...
        wls_m2_single(p_opt_m2, ValData, kinetic, nCutDOT_m2));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt_Modell2.mat'), 'p_opt_m2');

for e = 1:numel(TrainData)
    plot_fit_m2(TrainData{e}, p_opt_m2, kinetic, nCutDOT_m2, ...
                sprintf('Modell2 LHS | Training %s', TrainData{e}.Name));
end
plot_fit_m2(ValData, p_opt_m2, kinetic, nCutDOT_m2, ...
            sprintf('Modell2 LHS | Validierung %s', ValData.Name));


%% ======================================================================
%% LHS-Multistart
%% ======================================================================
function [p_opt, fval] = lhs_multistart(obj_fun, p0, pLB, pUB, options, N_lhs, K_opt)
% LHS-Multistart um fmincon. Startpunkte log-skaliert in [pLB,pUB].
% Voraussetzung: alle Grenzen > 0.
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
            [pk, Jk] = fmincon(obj_fun, Pstart(k,:), [], [], [], [], ...
                               pLBr, pUBr, [], options);
            if Jk < fval, fval = Jk;  p_opt = pk; end
        catch
        end
    end
    if wasCol, p_opt = p_opt(:); end
end


%% ======================================================================
%% WLS-Guetefunktionen Modell 1
%% ======================================================================
function J = wls_m1_multi(p, Train, kinetic, withOxygen, nCutDOT)
    J = 0;
    for e = 1:numel(Train)
        J = J + wls_m1_single(p, Train{e}, kinetic, withOxygen, nCutDOT);
    end
    if ~isfinite(J), J = 1e8; end
end

function J = wls_m1_single(p, D, kinetic, withOxygen, nCutDOT)
    M  = { D.Biomasse, 1, 'Biomasse'; D.Glucose, 2, 'Glucose' };
    x0 = [D.Biomasse.y(1); D.Glucose.y(1)];
    DOTstern = 0;
    if withOxygen
        DOTstern = max(D.O2.y);                 % aus vollem O2-Signal
        O2fit    = drop_leading(D.O2, nCutDOT); % fuehrende DOT-Punkte weg
        M        = [M; {O2fit, 3, 'DOT'}];
        x0       = [x0; O2fit.y(1)];            % Startwert = erster behaltener DOT
    end

    t_all = [];
    for i = 1:size(M,1), t_all = [t_all; M{i,1}.t(:)]; end
    t_all = unique(t_all);

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    try
        [~, X] = ode15s(@(t,x) Modell1(t,x,p,kinetic,withOxygen,DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:))), J = 1e8; return; end

    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1}; idxS = M{i,2}; name = M{i,3};
        [~, iT] = ismember(mess.t(:), t_all);
        r = (mess.y(:) - X(iT, idxS)) ./ sqrt(max(mess.var(:), eps));
        if name == "DOT", r = r / 20; end
        J = J + mean(r.^2);
    end
    if ~isfinite(J), J = 1e8; end
end


%% ======================================================================
%% WLS-Guetefunktionen Modell 2
%% ======================================================================
function J = wls_m2_multi(p, Train, kinetic, nCutDOT)
    J = 0;
    for e = 1:numel(Train)
        J = J + wls_m2_single(p, Train{e}, kinetic, nCutDOT);
    end
    if ~isfinite(J), J = 1e8; end
end

function J = wls_m2_single(p, D, kinetic, nCutDOT)
    O2fit = drop_leading(D.O2, nCutDOT);
    M = { D.Biomasse, 1, 'Biomasse'; ...
          D.Glucose,  2, 'Glucose';  ...
          D.Ammonium, 3, 'Ammonium'; ...
          D.Base,     4, 'Base';     ...
          O2fit,      5, 'DOT'       };
    x0 = [D.Biomasse.y(1); D.Glucose.y(1); D.Ammonium.y(1); D.Base.y(1); O2fit.y(1)];
    DOTstern = max(D.O2.y);

    t_all = [];
    for i = 1:size(M,1), t_all = [t_all; M{i,1}.t(:)]; end
    t_all = unique(t_all);

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    try
        [~, X] = ode15s(@(t,x) Modell2(t,x,p,kinetic,DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:))), J = 1e8; return; end

    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1}; idxS = M{i,2}; name = M{i,3};
        [~, iT] = ismember(mess.t(:), t_all);
        r = (mess.y(:) - X(iT, idxS)) ./ sqrt(max(mess.var(:), eps));
        if name == "DOT", r = r / 10; end
        J = J + mean(r.^2);
    end
    if ~isfinite(J), J = 1e8; end
end


%% ======================================================================
%% Hilfsfunktion: fuehrende Punkte entfernen
%% ======================================================================
function m = drop_leading(m, n)
% Entfernt die ersten n Punkte aus einer Messgroesse-Struct (Felder t,y,var).
    if n < 1, return; end
    keep = true(size(m.t));
    keep(1:min(n, numel(keep))) = false;
    m.t   = m.t(keep);
    m.y   = m.y(keep);
    m.var = m.var(keep);
end


%% ======================================================================
%% Visualisierung
%% ======================================================================
function plot_fit_m1(D, p, kinetic, withOxygen, nCutDOT, titel)
    x0 = [D.Biomasse.y(1); D.Glucose.y(1)];
    DOTstern = 0;  t_cutDOT = NaN;
    if withOxygen
        DOTstern = max(D.O2.y);
        O2fit    = drop_leading(D.O2, nCutDOT);
        x0       = [x0; O2fit.y(1)];
        if nCutDOT >= 1 && numel(D.O2.t) > nCutDOT
            t_cutDOT = D.O2.t(nCutDOT+1);
        end
    end

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
        plot(t_sim, X(:,3), 'LineWidth', 2);
        if ~isnan(t_cutDOT), xline(t_cutDOT, '--k', 'Cut', 'LabelVerticalAlignment','bottom'); end
        title('Sauerstoff (DOT)'); ylabel('cO_2 (%)');
        legend('Messung \pm \sigma','Simulation'); grid on;
    end
    xlabel('BatchAge (h)');
end


function plot_fit_m2(D, p, kinetic, nCutDOT, titel)
    O2fit = drop_leading(D.O2, nCutDOT);
    x0 = [D.Biomasse.y(1); D.Glucose.y(1); D.Ammonium.y(1); D.Base.y(1); O2fit.y(1)];
    DOTstern = max(D.O2.y);
    t_cutDOT = NaN;
    if nCutDOT >= 1 && numel(D.O2.t) > nCutDOT, t_cutDOT = D.O2.t(nCutDOT+1); end

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
        if k == 5 && ~isnan(t_cutDOT)
            xline(t_cutDOT, '--k', 'Cut', 'LabelVerticalAlignment','bottom');
        end
        title([titel ' - ' kanal{k,3}]); ylabel(kanal{k,4});
        legend('Messung \pm \sigma','Simulation'); grid on;
    end
    xlabel('BatchAge (h)');
end
