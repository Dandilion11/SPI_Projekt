%% main_04_parameter_fitting_LHS.m
%  Parameterfitting Modell 1 & 2 mit LHS-Multistart.
%  Training:    RamScDef03 (TrainData)
%  Validierung: RamScDef04 (ValData)
%  Fit ausschliesslich auf TrainData; ValData wird nur verglichen.
%
%  Features:
%   - LHS-Multistart um fmincon (billiges Screening -> beste Startpunkte).
%   - nCutDOT: fuehrende DOT-Punkte werden NICHT gefittet, x0(DOT) startet
%     beim ersten behaltenen Punkt (Anfangswert bleibt konsistent).
% =========================================================================
clear; clc; close all;

projectRoot = pwd;
load(fullfile(projectRoot, '..', 'Daten', 'Daten_Processed', 'Processed_Batch_Data.mat'));
addpath(fullfile(projectRoot, '..', 'Modelle'), '-begin');
rehash; clear Modell1 Modell2

kinetic    = 3;      % 3 = Monod
withOxygen = true;

% LHS-Multistart-Konfiguration
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
N_lhs = 100;   % LHS-Screening-Punkte (zum Testen z.B. 30)
K_opt = 3;     % beste Startpunkte -> teurer fmincon (zum Testen z.B. 2)

%% ======================================================================
%% MODELL 1  (Training auf RamScDef03)
%% ======================================================================
nCutDOT = 1;    % fuehrende DOT-Punkte, die NICHT gefittet werden

p0  = [0.3;  0.5;  0.15];
pLB = [0.01; 0.01; 0.01];
pUB = [1.0;  500;  1.0];
if withOxygen
    p0  = [p0;  1.0; 200];
    pLB = [pLB; 0.01;   1];
    pUB = [pUB; 500; 1000];
end

[TrainFit1, x0_tr1] = prep_m1(TrainData, nCutDOT, withOxygen);
obj_m1 = @(p) calculate_wls_error(p, x0_tr1, TrainFit1, kinetic, withOxygen);
[p_opt, fval] = lhs_multistart(obj_m1, p0, pLB, pUB, options, N_lhs, K_opt);

fprintf('\n--- Modell1 (LHS): Fitting abgeschlossen (Training RamScDef03) ---\n');
fprintf('WLS-Fehler (Training):    %.4f\n', fval);
fprintf('mu_max = %.4f 1/h\n', p_opt(1));
fprintf('K_S    = %.4f g/L\n', p_opt(2));
fprintf('Y_XS   = %.4f g/g\n', p_opt(3));
if withOxygen
    fprintf('Y_XO   = %.4f g/g\n', p_opt(4));
    fprintf('KLa    = %.4f 1/h\n', p_opt(5));
end
[ValFit1, x0_val1] = prep_m1(ValData, nCutDOT, withOxygen);
fprintf('WLS-Fehler (Validierung RamScDef04): %.4f\n', ...
        calculate_wls_error(p_opt, x0_val1, ValFit1, kinetic, withOxygen));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt.mat'), 'p_opt');

plot_m1(TrainData, p_opt, kinetic, withOxygen, nCutDOT, 'Modell1 LHS | Training RamScDef03');
plot_m1(ValData,   p_opt, kinetic, withOxygen, nCutDOT, 'Modell1 LHS | Validierung RamScDef04');


%% ======================================================================
%% MODELL 2  (Training auf RamScDef03)
%% ======================================================================
nCutDOT_m2 = 1;

p0_m2  = [0.3;  0.5;  0.15; 1.0;  0.05;  1.0; 200];
pLB_m2 = [0.01; 0.01; 0.01; 0.01; 0.001; 0.01;   1];
pUB_m2 = [1.0;  5.0;  1.0;  10;   1.0;   100; 1000];

[TrainFit2, x0_tr2] = prep_m2(TrainData, nCutDOT_m2);
obj_m2 = @(p) calculate_wls_error_m2(p, x0_tr2, TrainFit2, kinetic);
[p_opt_m2, fval_m2] = lhs_multistart(obj_m2, p0_m2, pLB_m2, pUB_m2, options, N_lhs, K_opt);

fprintf('\n--- Modell2 (LHS): Fitting abgeschlossen (Training RamScDef03) ---\n');
fprintf('WLS-Fehler (Training):    %.4f\n', fval_m2);
fprintf('mu_max = %.4f 1/h\n',  p_opt_m2(1));
fprintf('K_S    = %.4f g/L\n',  p_opt_m2(2));
fprintf('Y_XS   = %.4f g/g\n',  p_opt_m2(3));
fprintf('Y_Bam  = %.4f\n',      p_opt_m2(4));
fprintf('Y_AmX  = %.4f\n',      p_opt_m2(5));
fprintf('Y_XO   = %.4f g/g\n',  p_opt_m2(6));
fprintf('KLa    = %.4f 1/h\n',  p_opt_m2(7));
[ValFit2, x0_val2] = prep_m2(ValData, nCutDOT_m2);
fprintf('WLS-Fehler (Validierung RamScDef04): %.4f\n', ...
        calculate_wls_error_m2(p_opt_m2, x0_val2, ValFit2, kinetic));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt_Modell2.mat'), 'p_opt_m2');

plot_m2(TrainData, p_opt_m2, kinetic, nCutDOT_m2, 'Modell2 LHS | Training RamScDef03');
plot_m2(ValData,   p_opt_m2, kinetic, nCutDOT_m2, 'Modell2 LHS | Validierung RamScDef04');


%% ======================================================================
%% Vorbereitung (nCutDOT anwenden, x0 bauen)
%% ======================================================================
function [D, x0] = prep_m1(Data, nCutDOT, withOxygen)
    D  = Data;
    x0 = [Data.Biomasse.y(1); Data.Glucose.y(1)];
    if withOxygen
        x0 = [x0; Data.O2.y(1)];
        if nCutDOT >= 1 && numel(Data.O2.t) > nCutDOT
            D.O2    = drop_leading(Data.O2, nCutDOT);
            x0(end) = Data.O2.y(nCutDOT+1);
        end
    end
end

function [D, x0] = prep_m2(Data, nCutDOT)
    D  = Data;
    x0 = [Data.Biomasse.y(1); Data.Glucose.y(1); Data.Ammonium.y(1); Data.Base.y(1); Data.O2.y(1)];
    if nCutDOT >= 1 && numel(Data.O2.t) > nCutDOT
        D.O2    = drop_leading(Data.O2, nCutDOT);
        x0(end) = Data.O2.y(nCutDOT+1);
    end
end

function m = drop_leading(m, n)
% Entfernt die ersten n Punkte einer Messgroesse (Felder t,y,var).
    if n < 1, return; end
    keep = true(size(m.t));
    keep(1:min(n, numel(keep))) = false;
    m.t = m.t(keep);  m.y = m.y(keep);  m.var = m.var(keep);
end


%% ======================================================================
%% LHS-Multistart
%% ======================================================================
function [p_opt, fval] = lhs_multistart(obj_fun, p0, pLB, pUB, options, N_lhs, K_opt)
% Log-skalierte LHS-Startpunkte in [pLB,pUB], billiges Screening, dann
% fmincon von p0 + den K_opt besten. Voraussetzung: alle Grenzen > 0.
    wasCol = iscolumn(p0);
    p0r = p0(:).';  pLBr = pLB(:).';  pUBr = pUB(:).';
    d   = numel(p0r);

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
%% Visualisierung
%% ======================================================================
function plot_m1(Data, p, kinetic, withOxygen, nCutDOT, titel)
    x0 = [Data.Biomasse.y(1); Data.Glucose.y(1)];
    DOTstern = 0;  t_cutDOT = NaN;
    if withOxygen
        DOTstern = max(Data.O2.y);
        x0 = [x0; Data.O2.y(1)];
        if nCutDOT >= 1 && numel(Data.O2.t) > nCutDOT
            x0(end)  = Data.O2.y(nCutDOT+1);
            t_cutDOT = Data.O2.t(nCutDOT+1);
        end
    end

    t_end = max([Data.Biomasse.t(:); Data.Glucose.t(:)]) + 1;
    t_sim = linspace(0, t_end, 200);
    [~, X] = ode15s(@(t,x) Modell1(t,x,p,kinetic,withOxygen,DOTstern), ...
                    t_sim, x0, odeset('RelTol',1e-5,'AbsTol',1e-7));

    n = 2 + double(withOxygen);
    figure('Name', titel, 'Position', [200 150 850 600]);
    subplot(n,1,1);
    errorbar(Data.Biomasse.t, Data.Biomasse.y, sqrt(Data.Biomasse.var), 'o', 'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, X(:,1), 'LineWidth', 2); title([titel ' - Biomasse']);
    ylabel('c_X (g/L)'); legend('Messung \pm \sigma','Simulation'); grid on;
    subplot(n,1,2);
    errorbar(Data.Glucose.t, Data.Glucose.y, sqrt(Data.Glucose.var), 'o', 'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, X(:,2), 'LineWidth', 2); title('Glucose');
    ylabel('c_{Glc} (g/L)'); legend('Messung \pm \sigma','Simulation'); grid on;
    if withOxygen
        subplot(n,1,3);
        errorbar(Data.O2.t, Data.O2.y, sqrt(Data.O2.var), 'o', 'MarkerFaceColor','b','MarkerSize',4); hold on;
        plot(t_sim, X(:,3), 'LineWidth', 2);
        if ~isnan(t_cutDOT), xline(t_cutDOT, '--k', 'Cut', 'LabelVerticalAlignment','bottom'); end
        title('Sauerstoff (DOT)'); ylabel('cO_2 (%)');
        legend('Messung \pm \sigma','Simulation'); grid on;
    end
    xlabel('BatchAge (h)');
end

function plot_m2(Data, p, kinetic, nCutDOT, titel)
    x0 = [Data.Biomasse.y(1); Data.Glucose.y(1); Data.Ammonium.y(1); Data.Base.y(1); Data.O2.y(1)];
    DOTstern = max(Data.O2.y);
    t_cutDOT = NaN;
    if nCutDOT >= 1 && numel(Data.O2.t) > nCutDOT
        x0(end)  = Data.O2.y(nCutDOT+1);
        t_cutDOT = Data.O2.t(nCutDOT+1);
    end

    t_end = max([Data.Biomasse.t(:); Data.Glucose.t(:); Data.Ammonium.t(:); Data.Base.t(:); Data.O2.t(:)]) + 1;
    t_sim = linspace(0, t_end, 300);
    [~, X] = ode15s(@(t,x) Modell2(t,x,p,kinetic,DOTstern), ...
                    t_sim, x0, odeset('RelTol',1e-5,'AbsTol',1e-7));

    kanal = {Data.Biomasse,1,'Biomasse','c_X (g/L)'; ...
             Data.Glucose, 2,'Glucose', 'c_{Glc} (g/L)'; ...
             Data.Ammonium,3,'Ammonium','c_{Am} (g/L)'; ...
             Data.Base,    4,'Base',    'c_{Base}'; ...
             Data.O2,      5,'DOT',     'cO_2 (%)'};
    figure('Name', titel, 'Position', [250 60 900 900]);
    for k = 1:5
        subplot(5,1,k);
        errorbar(kanal{k,1}.t, kanal{k,1}.y, sqrt(kanal{k,1}.var), 'o', 'MarkerFaceColor','b','MarkerSize',4); hold on;
        plot(t_sim, X(:, kanal{k,2}), 'LineWidth', 2);
        if k == 5 && ~isnan(t_cutDOT), xline(t_cutDOT, '--k', 'Cut', 'LabelVerticalAlignment','bottom'); end
        title([titel ' - ' kanal{k,3}]); ylabel(kanal{k,4});
        legend('Messung \pm \sigma','Simulation'); grid on;
    end
    xlabel('BatchAge (h)');
end


%% ======================================================================
%% WLS-Guetefunktionen
%% ======================================================================
function J = calculate_wls_error(p, x0, Data, kinetic, withOxygen)
% WLS Modell 1. Zustaende ohne O2: [cX; cGlc], mit O2: [cX; cGlc; DOT].
    M = { Data.Biomasse, 1, 'Biomasse'; Data.Glucose, 2, 'Glucose' };
    DOTstern = 0;
    if withOxygen
        M = [M; {Data.O2, 3, 'DOT'}];
        DOTstern = max(Data.O2.y);
    end
    wmode = 'mean';
    wsig  = ones(size(M,1),1);

    t_all = [];
    for i = 1:size(M,1), t_all = [t_all; M{i,1}.t(:)]; end
    t_all = unique(t_all);

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    try
        [~, X_sim] = ode15s(@(t,x) Modell1(t,x,p,kinetic,withOxygen,DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X_sim,1) ~= numel(t_all) || any(~isfinite(X_sim(:))), J = 1e8; return; end

    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1}; idxState = M{i,2}; name = M{i,3};
        [tf, iT] = ismember(mess.t(:), t_all);
        if any(~tf), warning('Messzeitpunkt fuer %s nicht gefunden.', name); J = 1e8; return; end
        y_sim = X_sim(iT, idxState);
        var_i = max(mess.var(:), eps);
        r = (mess.y(:) - y_sim) ./ sqrt(var_i);
        if name == "DOT", r = r/20; end
        switch wmode
            case 'sum',  contrib = sum(r.^2);
            case 'mean', contrib = mean(r.^2);
        end
        J = J + wsig(i) * contrib;
    end
    if ~isfinite(J), J = 1e8; end
end

function J = calculate_wls_error_m2(p, x0, Data, kinetic)
% WLS Modell 2. Zustaende: [cX; cGlc; cAm; cBase; DOT].
    M = { Data.Biomasse, 1, 'Biomasse'; ...
          Data.Glucose,  2, 'Glucose';  ...
          Data.Ammonium, 3, 'Ammonium'; ...
          Data.Base,     4, 'Base';     ...
          Data.O2,       5, 'DOT'       };
    DOTstern = max(Data.O2.y);
    wmode = 'mean';
    wsig  = [1; 1; 1; 1; 1];

    t_all = [];
    for i = 1:size(M,1), t_all = [t_all; M{i,1}.t(:)]; end
    t_all = unique(t_all);

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    try
        [~, X_sim] = ode15s(@(t,x) Modell2(t,x,p,kinetic,DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X_sim,1) ~= numel(t_all) || any(~isfinite(X_sim(:))), J = 1e8; return; end

    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1}; idxState = M{i,2}; name = M{i,3};
        [tf, iT] = ismember(mess.t(:), t_all);
        if any(~tf), warning('Messzeitpunkt fuer %s nicht gefunden.', name); J = 1e8; return; end
        y_sim = X_sim(iT, idxState);
        var_i = max(mess.var(:), eps);
        r = (mess.y(:) - y_sim) ./ sqrt(var_i);
        if name == "DOT", r = r/10; end
        switch wmode
            case 'sum',  contrib = sum(r.^2);
            case 'mean', contrib = mean(r.^2);
        end
        J = J + wsig(i) * contrib;
    end
    if ~isfinite(J), J = 1e8; end
end
