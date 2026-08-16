% main_04_parameter_fitting_m3_woEtOH_LHS_10p_multi.m
%
% Parameteridentifikation Modell 3 (Fed-Batch, ohne Ethanol).
% EIN gemeinsamer Parametersatz fuer alle vier Trainingsexperimente,
% Validierung auf einem zurueckgehaltenen Lauf.
%
%   - Mittelung pro Messkanal statt Aufsummieren
%   - Base wird als Absolutwert bewertet (Kraemer & King 2017, Gl. 22)
%   - KLa, KAm, KPh sind fixiert -> 7 freie Parameter
%   - Optimierung in log-Skala, Multistart per LHS

clear; clc; close all;
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

%% 1. Daten laden ---------------------------------------------------------
projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fullfile(projectRoot,'..', '..');
load(fullfile(projectRoot,'Daten/Daten_Processed/Processed_FedBatch_Modell3_MultiExp.mat'));
addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),'-begin');
rehash;

SIMFUN = @(tt,xx,uu,pp,dd,PP) sim_m3_sample_10p(tt,xx,uu,pp,dd,PP,'fast');

nTrain = numel(TrainSet);
nVal   = numel(ValSet);

DOTstern_train = arrayfun(@(D) max(D.O2.y), TrainSet);
DOTstern_val   = arrayfun(@(D) max(D.O2.y), ValSet);

fprintf('Training   : %s\n', strjoin({TrainSet.name}, ', '));
fprintf('Validierung: %s\n', strjoin({ValSet.name}, ', '));

%% 2. Kanalgewichte -------------------------------------------------------
% [Biomasse Glucose Ammonium Phosphat Base DOT], 0 schaltet einen Kanal ab.
wsig = [1 1 1 1 0.01 0.0001];

%% 3. Parameter -----------------------------------------------------------
namen = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO','KAm','KPh'};
p0  = [0.328 4.32 0.1412 0.087 0.087 0.0245 500 0.47 0.0001 0.0001];
pLB = [ 0.01,  1e-3,  0.05,  1e-3,  1e-4,  1e-4,  50,   0.01,  0.0001,  0.0001];
pUB = [ 1.00,  10.0,  1.00,  1.00,  1.00,  1.0,   800,  10,    1.0,     0.5   ];

% Fixierte Parameter -- pLB == pUB haelt sie in fmincon fest.
%   KLa      : nur das Produkt KLa*YXO ist bestimmbar (DOT quasistationaer).
%   KAm, KPh : nicht identifizierbar, bleiben aber im Modell -- ohne sie
%              laufen mAm/mPh negativ. Werte nach Herold et al. 2017.
iFix = [7 9 10];
pFix = [500, 0.0001, 0.0001];
p0(iFix) = pFix;  pLB(iFix) = pFix;  pUB(iFix) = pFix;
isFree = true(1,numel(p0));  isFree(iFix) = false;

fprintf('Frei: %s\n', strjoin(namen(isFree), ', '));
fprintf('Fix : %s\n', strjoin(compose('%s=%.4g', string(namen(iFix))', pFix'), ', '));

% Bisher bestes Ergebnis: Startpunkt UND Akzeptanzschwelle.
pRef = [0.328 4.32 0.1412 0.087 0.087 0.0245 500 0.47 0.0001 0.0001];

options = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp', ...
                       'MaxFunctionEvaluations', 8000, ...
                       'FiniteDifferenceType','central', ...
                       'FiniteDifferenceStepSize', 1e-6, ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance', 1e-10);

obj_fun = @(p) wls_error_multi(p, TrainSet, DOTstern_train, SIMFUN, wsig);

% Optimiert wird ueber q = log10(p): die Parameter spannen vier
% Groessenordnungen, in log-Skala ist ein Schritt fuer alle dieselbe
% RELATIVE Aenderung. Das Optimum bleibt unveraendert.
obj_log = @(q) obj_fun(10.^q);
qLB = log10(pLB);   qUB = log10(pUB);

J_ref = obj_fun(pRef);
fprintf('\nReferenzpunkt pRef: J = %.4f  (Endergebnis muss <= sein)\n', J_ref);

%% 4. Startpunkte per Latin Hypercube Sampling ---------------------------
N_lhs = 30;
K_opt = 1;

rng(42);
L = lhsdesign(N_lhs, numel(p0));
P = 10.^(log10(pLB) + L .* (log10(pUB) - log10(pLB)));   % log-gleichverteilt
P(:,iFix) = repmat(pFix, N_lhs, 1);

fprintf('\nLHS-Screening ueber %d Punkte ...\n', N_lhs);
Jscreen = inf(N_lhs,1);
for k = 1:N_lhs
    try, Jscreen(k) = obj_fun(P(k,:)); catch, end
end
[~, order] = sort(Jscreen);
nGood = min(K_opt, sum(isfinite(Jscreen)));

Pstart = unique([pRef; p0; P(order(1:nGood), :)], 'rows', 'stable');
Qstart = log10(Pstart);

fprintf('Bester LHS-Punkt: J = %.4f\n', min(Jscreen));

%% 5. Optimierung ---------------------------------------------------------
p_opt = pRef;  fval = J_ref;
for k = 1:size(Qstart,1)
    fprintf('\n--- fmincon Start %d/%d ---\n', k, size(Qstart,1));
    try
        [qk, Jk] = fmincon(obj_log, Qstart(k,:), [], [], [], [], ...
                           qLB, qUB, [], options);
        fprintf('Start %d: J = %.4f\n', k, Jk);
        if Jk < fval, fval = Jk;  p_opt = 10.^qk; end
    catch Me
        fprintf('Start %d fehlgeschlagen: %s\n', k, Me.message);
    end
end

%% 6. Ergebnis ------------------------------------------------------------
fprintf('\n=====================================================\n');
fprintf('  Modell3 woEtOH -- Multi-Experiment-Fit\n');
fprintf('  %d freie Parameter | wsig = %s\n', nnz(isFree), mat2str(wsig));
fprintf('=====================================================\n');
fprintf('WLS gesamt (Training): %.4f   (Referenz war %.4f)\n\n', fval, J_ref);

% Klebt ein freier Parameter an einer Grenze, ist das Ergebnis von der
% Grenze bestimmt und die FIM dort nicht aussagekraeftig.
onBound = isFree & (p_opt <= pLB*1.001 | p_opt >= pUB*0.999);
for i = 1:numel(p_opt)
    if ~isFree(i)
        fprintf('%-8s = %10.4f   (fixiert)\n', namen{i}, p_opt(i));
    elseif onBound(i)
        fprintf('%-8s = %10.4f  <-- AN GRENZE\n', namen{i}, p_opt(i));
    else
        fprintf('%-8s = %10.4f\n', namen{i}, p_opt(i));
    end
end
fprintf('  -> Y_XO ist Y_XO_eff = KLa*YXO (KLa fixiert auf %.1f)\n', p_opt(7));

fprintf('\n--- Beitrag pro Experiment ---\n');
Jtr = zeros(nTrain,1);
for k = 1:nTrain
    [Jtr(k), contrib, nch] = wls_error_m3(p_opt, TrainSet(k), ...
                             DOTstern_train(k), SIMFUN, wsig);
    fprintf('%-12s J = %10.4f | chi2/N_eff = %6.1f\n', ...
            TrainSet(k).name, Jtr(k), Jtr(k)/max(nch,1));
    print_channels(contrib);
end
fprintf('\nchi2/N_eff gesamt (Training) = %.1f   (Erwartung ~1)\n', ...
        fval / (nnz(wsig)*nTrain));

%% 7. Lokale Sensitivitaet ------------------------------------------------
% Wie stark reagiert J auf +-20 % je Parameter? Kleine Werte heissen:
% mit DIESEN Daten nicht bestimmbar -- das motiviert die Versuchsplanung.
fprintf('\n--- Lokale Sensitivitaet ---\n');
fprintf('%-8s %12s %12s %12s\n', 'Param', 'J(-20%)', 'J(+20%)', 'max dJ [%]');
for j = 1:numel(p_opt)
    if ~isFree(j)
        fprintf('%-8s %12s %12s %12s   (fixiert)\n', namen{j}, '-','-','-');
        continue
    end
    pm = p_opt; pm(j) = max(p_opt(j)*0.8, pLB(j));
    pp = p_opt; pp(j) = min(p_opt(j)*1.2, pUB(j));
    Jm = obj_fun(pm);  Jp = obj_fun(pp);
    rel = 100 * max(abs([Jm Jp] - fval)) / fval;
    marker = ''; if rel < 1, marker = '   <-- unempfindlich'; end
    fprintf('%-8s %12.2f %12.2f %12.2f%s\n', namen{j}, Jm, Jp, rel, marker);
end

%% 8. Validierung ---------------------------------------------------------
fprintf('\n--- Validierung ---\n');
Jval = zeros(nVal,1);
for k = 1:nVal
    [Jval(k), contrib, nch] = wls_error_m3(p_opt, ValSet(k), ...
                              DOTstern_val(k), SIMFUN, wsig);
    fprintf('%-12s J = %10.4f | chi2/N_eff = %6.1f\n', ...
            ValSet(k).name, Jval(k), Jval(k)/max(nch,1));
    print_channels(contrib);
end
fprintf('\nVerhaeltnis Val/Train = %.2f\n', ...
        (sum(Jval)/nVal) / (sum(Jtr)/nTrain));

%% 9. Speichern -----------------------------------------------------------
saveDir = fullfile(projectRoot,'Daten','p_opt');
if ~exist(saveDir,'dir'); mkdir(saveDir); end
trainNames = {TrainSet.name};  valNames = {ValSet.name};
save(fullfile(saveDir,'p_opt_Modell3_woEtOH_10p_multi.mat'), ...
     'p_opt','namen','pLB','pUB','fval','Jtr','Jval','wsig','iFix','pFix', ...
     'isFree','J_ref','trainNames','valNames','DOTstern_train','DOTstern_val');
fprintf('\nGespeichert: %s\n', saveDir);

%% 10. Plots --------------------------------------------------------------
for k = 1:nTrain
    plot_experiment(TrainSet(k), p_opt, DOTstern_train(k), SIMFUN, ...
                    sprintf('Training: %s', TrainSet(k).name));
end
for k = 1:nVal
    plot_experiment(ValSet(k), p_opt, DOTstern_val(k), SIMFUN, ...
                    sprintf('Validierung: %s', ValSet(k).name));
end

bildordner = fullfile(projectRoot,'Bilder','Modell3');
figs = findobj('Type','figure');
[~, ord] = sort([figs.Number]);
for k = ord(:).'
    save_fig(figs(k), figs(k).Name, bildordner);
end


%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function J = wls_error_multi(p, Set, DOTstern, SIMFUN, wsig)
% Summe des Guetefunktionals ueber alle Experimente.
    J = 0;
    for k = 1:numel(Set)
        J = J + wls_error_m3(p, Set(k), DOTstern(k), SIMFUN, wsig);
        if ~isfinite(J), J = 1e8; return; end
    end
end


function [J, contrib, nch] = wls_error_m3(p, D, DOTstern, SIMFUN, wsig)
% Guetefunktional fuer EIN Experiment.
%   Zustaende:  x = [V; mX; mGlc; mAm; mPh; mB; DOT]
%   Residuum:   r = (y - y_sim)/sigma
%   Gewichtung: Mittelwert pro Kanal, sonst dominieren Base und DOT mit
%               ihren hunderten Punkten. Dadurch ist N_eff = Anzahl
%               aktiver Kanaele.

    % {Messreihe, Zustandsindex, durch V teilen?}
    M = { D.Biomasse, 2, true;   D.Glucose,  3, true; ...
          D.Ammonium, 4, true;   D.Phosphat, 5, true; ...
          D.Base,     6, false;  D.O2,       7, false };

    contrib = zeros(1,6);  nch = 0;

    t_all = [];
    for i = 1:6
        if wsig(i) > 0, t_all = [t_all; M{i,1}.t(:)]; end %#ok<AGROW>
    end
    if isempty(t_all), J = 1e8; return; end
    tu = D.u(1,:).';
    t_all = unique([t_all; tu(tu >= min(t_all) & tu <= max(t_all))]);

    try
        X = SIMFUN(t_all, D.x0, D.u, p, DOTstern, D.Probe);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:)))
        J = 1e8; return;
    end

    V = X(:,1);  J = 0;
    for i = 1:6
        if wsig(i) == 0, continue; end
        mess = M{i,1};
        if isempty(mess.t), continue; end

        [~, iT] = ismember(mess.t(:), t_all);
        y_sim = X(iT, M{i,2});
        if M{i,3}, y_sim = y_sim ./ V(iT); end

        r = (mess.y(:) - y_sim) ./ sqrt(max(mess.var(:), eps));

        contrib(i) = wsig(i) * mean(r.^2);
        J          = J + contrib(i);
        nch        = nch + 1;
    end
end


function print_channels(contrib)
    nm = {'Biomasse','Glucose','Ammonium','Phosphat','Base','DOT'};
    fprintf('     ');
    for i = 1:numel(nm)
        fprintf('%s=%.1f  ', nm{i}, contrib(i));
    end
    fprintf('\n');
end


function plot_experiment(D, p_opt, DOTstern, SIMFUN, titel)
% Simulation ueber den vollen Horizont, alle sechs Kanaele uebereinander.
    FS = 14;  FSA = 14;

    t_sim = linspace(D.u(1,1), max([D.Biomasse.t; D.O2.t; D.Base.t])+1, 400);
    X = SIMFUN(t_sim, D.x0, D.u, p_opt, DOTstern, D.Probe);
    V = X(:,1);

    figure('Name', titel, 'Position', [200 40 950 1000]);
    plot_row(1, D.Biomasse, t_sim, X(:,2)./V, 'Biomasse', 'c_X (g/L)', FS, FSA);
    plot_row(2, D.Glucose,  t_sim, X(:,3)./V, 'Glucose',  'c_{Glc} (g/L)', FS, FSA);
    plot_row(3, D.Ammonium, t_sim, X(:,4)./V, 'Ammonium', 'c_{Am} (g/L)', FS, FSA);
    plot_row(4, D.Phosphat, t_sim, X(:,5)./V, 'Phosphat', 'c_{Ph} (g/L)', FS, FSA);
    plot_row(5, D.Base,     t_sim, X(:,6),    'Base',     'm_B (L)', FS, FSA);
    plot_row(6, D.O2,       t_sim, X(:,7),    'DOT',      'DOT (%)', FS, FSA);
    xlabel('BatchAge (h)', 'FontSize', FS);
    sgtitle(titel, 'FontSize', FS);
end


function plot_row(row, mess, t_sim, y_sim, name, ylab, FS, FSA)
    subplot(6,1,row);
    errorbar(mess.t, mess.y, sqrt(mess.var), 'o', ...
             'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, y_sim, 'LineWidth', 2);
    title(name, 'FontSize', FS); ylabel(ylab, 'FontSize', FS);
    legend('Messung ± σ','Simulation','Location','best', 'FontSize', FS);
    set(gca, 'FontSize', FSA);
    grid on;
end


function save_fig(fig, name, ordner)
    if ~exist(ordner,'dir'), mkdir(ordner); end
    name = regexprep(strtrim(name), '[^\w\-]', '_');
    if isempty(name), name = sprintf('Figure_%d', fig.Number); end
    basis = fullfile(ordner, name);

    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
    drawnow;

    try
        exportgraphics(fig, [basis '.svg'], 'ContentType', 'vector');
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, [basis '.svg'], '-dsvg');
    end

    fprintf('  gespeichert: %s.svg\n', basis);
end