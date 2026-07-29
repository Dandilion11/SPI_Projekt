%% main_04_parameter_fitting.m
%  Parameterfitting fuer Modell 1 & 2 auf MEHREREN Trainingsexperimenten
%  (ein gemeinsamer Parametersatz) + Validierung auf einem Experiment.
%  Nicht-LHS-Variante, Basisausgang.
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

%% ======================================================================
%% MODELL 1
%% ======================================================================
withOxygen = true;

% Parametervektor [mu_max; K_S; Y_XS] (+ [Y_XO; KLa] bei withOxygen)
p0  = [0.3;  0.5;  0.15];
pLB = [0.01; 0.01; 0.01];
pUB = [1.0;  500;  1.0];
if withOxygen
    p0  = [p0;  1.0; 200];
    pLB = [pLB; 0.01;  1];
    pUB = [pUB; 500; 1000];
end

obj_m1  = @(p) wls_m1_multi(p, TrainData, kinetic, withOxygen);
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[p_opt, fval] = fmincon(obj_m1, p0, [], [], [], [], pLB, pUB, [], options);

fprintf('\n--- Modell1: Parameteridentifikation abgeschlossen ---\n');
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
        wls_m1_single(p_opt, ValData, kinetic, withOxygen));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt.mat'), 'p_opt');

% Basis-Visualisierung: jedes Trainingsexperiment + Validierung
for e = 1:numel(TrainData)
    plot_fit_m1(TrainData{e}, p_opt, kinetic, withOxygen, ...
                sprintf('Modell1 | Training %s', TrainData{e}.Name));
end
plot_fit_m1(ValData, p_opt, kinetic, withOxygen, ...
            sprintf('Modell1 | Validierung %s', ValData.Name));


%% ======================================================================
%% MODELL 2
%% ======================================================================
% Parameter [mu_max; K_S; Y_XS; Y_Bam; Y_AmX; Y_XO; KLa]
p0_m2  = [0.3;  0.5;  0.15; 1.0;  0.05;  1.0; 200];
pLB_m2 = [0.01; 0.01; 0.01; 0.01; 0.001; 0.01;   1];
pUB_m2 = [1.0;  5.0;  1.0;  10;   1.0;   100; 1000];

obj_m2 = @(p) wls_m2_multi(p, TrainData, kinetic);
[p_opt_m2, fval_m2] = fmincon(obj_m2, p0_m2, [], [], [], [], pLB_m2, pUB_m2, [], options);

fprintf('\n--- Modell2: Parameteridentifikation abgeschlossen ---\n');
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
        wls_m2_single(p_opt_m2, ValData, kinetic));

save(fullfile(projectRoot, '..', 'Daten', 'p_opt_Modell2.mat'), 'p_opt_m2');

for e = 1:numel(TrainData)
    plot_fit_m2(TrainData{e}, p_opt_m2, kinetic, ...
                sprintf('Modell2 | Training %s', TrainData{e}.Name));
end
plot_fit_m2(ValData, p_opt_m2, kinetic, ...
            sprintf('Modell2 | Validierung %s', ValData.Name));


%% ======================================================================
%% WLS-Guetefunktionen Modell 1
%% ======================================================================
function J = wls_m1_multi(p, Train, kinetic, withOxygen)
% Summe der WLS-Fehler ueber alle Trainingsexperimente (gemeinsames p).
    J = 0;
    for e = 1:numel(Train)
        J = J + wls_m1_single(p, Train{e}, kinetic, withOxygen);
    end
    if ~isfinite(J), J = 1e8; end
end

function J = wls_m1_single(p, D, kinetic, withOxygen)
% WLS eines einzelnen Experiments fuer Modell 1.
% Zustaende ohne O2: [cX; cGlc], mit O2: [cX; cGlc; DOT].
% x0 aus der jeweils ersten Messung dieses Experiments.
    M  = { D.Biomasse, 1, 'Biomasse'; D.Glucose, 2, 'Glucose' };
    x0 = [D.Biomasse.y(1); D.Glucose.y(1)];
    DOTstern = 0;
    if withOxygen
        M  = [M; {D.O2, 3, 'DOT'}];
        x0 = [x0; D.O2.y(1)];
        DOTstern = max(D.O2.y);
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
        if name == "DOT", r = r / 20; end   % viele DOT-Punkte abdaempfen
        J = J + mean(r.^2);                 % jede Messgroesse gleich gewichtet
    end
    if ~isfinite(J), J = 1e8; end
end


%% ======================================================================
%% WLS-Guetefunktionen Modell 2
%% ======================================================================
function J = wls_m2_multi(p, Train, kinetic)
    J = 0;
    for e = 1:numel(Train)
        J = J + wls_m2_single(p, Train{e}, kinetic);
    end
    if ~isfinite(J), J = 1e8; end
end

function J = wls_m2_single(p, D, kinetic)
% WLS eines einzelnen Experiments fuer Modell 2.
% Zustaende: [cX; cGlc; cAm; cBase; DOT].
    M = { D.Biomasse, 1, 'Biomasse'; ...
          D.Glucose,  2, 'Glucose';  ...
          D.Ammonium, 3, 'Ammonium'; ...
          D.Base,     4, 'Base';     ...
          D.O2,       5, 'DOT'       };
    x0 = [D.Biomasse.y(1); D.Glucose.y(1); D.Ammonium.y(1); D.Base.y(1); D.O2.y(1)];
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
        if name == "DOT", r = r / 20; end
        J = J + mean(r.^2);
    end
    if ~isfinite(J), J = 1e8; end
end


%% ======================================================================
%% Visualisierung
%% ======================================================================
function plot_fit_m1(D, p, kinetic, withOxygen, titel)
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


function plot_fit_m2(D, p, kinetic, titel)
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
