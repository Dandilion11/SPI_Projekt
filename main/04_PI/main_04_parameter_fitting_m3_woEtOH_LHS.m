% main_04d_parameter_fitting_m3_woEtOH_fromfeed.m
% Wie main_04 (Modell3_woEtOH), aber gefittet wird erst AB der ersten
% Fuetterung. Die vielen dichten Baseline-Punkte vor dem ersten Feed
% (flach, ~0, kaum Info ueber die Kinetik) werden aus der Zielfunktion
% entfernt, damit sie den Fit nicht ueberproportional dominieren.
%
% Wichtig: die ODE integriert weiterhin ab x0 bei t_start = u(1,1).
% Es werden NUR die Messpunkte mit t < t_feed nicht mehr bewertet -
% die Anfangsbedingung bleibt also korrekt.
clear; clc; close all;

% 1. Load Preprocessed Data
projectRoot = pwd;
load(fullfile(projectRoot,'..','Daten/Daten_Processed/Processed_FedBatch_Modell3.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
addpath(fullfile(projectRoot, '..','utils'),'-begin');
rehash; clear Modell3_woEtOH

Data = TrainData;
Probe_train = TrainDataProbe;
u    = Data.u;
x0   = Data.x0;                % [V; mX; mGlc; mAm; mPh; mB; DOT]

% ---- Ersten Feed-Zeitpunkt bestimmen ----------------------------------
% Substrat-Feedraten in u: Zeile 2 = uAm, 4 = uPh, 6 = uGlc.
% Erster Zeitpunkt, zu dem irgendeine davon > 0 wird.
feedRows   = [2 4 6];
feedOn     = any(u(feedRows,:) > 0, 1);
idxFeed    = find(feedOn, 1, 'first');
if isempty(idxFeed)
    t_feed = u(1,1);                 % Fallback: kein Feed gefunden
else
    t_feed = u(1, idxFeed);
end


% t_feed = 0.0;   % <- bei Bedarf manuell ueberschreiben
t_hi = 60;        % 60 / bei inf einsetzen für hinteres kappen ausschalten
DataFull = Data;                             % voller Satz nur fuer den Plot
DataFit  = window_data(Data, t_feed, t_hi);   % inf / nur Punkte ab t_feed fuers Fitting
fprintf('Erster Feed bei t_feed = %.3f h. Fitting nur fuer t >= t_feed.\n', t_feed);
report_window(DataFull, DataFit);
% -----------------------------------------------------------------------


% 2. Anfangswerte und Parameter
namen = {'mumax','KS','YXS', 'YAmX','YPhX','YB_Am','KLa','YXO'};
p0 =    [0.30,   0.50, 0.15,  0.05,  0.02,  1.0,     200,  1.0];
pLB =   [0.01,   0.01, 0.05, 0.001, 0.001,  0.1,     10,   0.1];
pUB =   [1.00,   5.00, 1.00, 1.000, 1.000,  5.0,     800,  5.0];

% Parameteridentifikation (WLS) -- auf DataFit (ab Feed)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', 5000, ...
                       'FiniteDifferenceType', 'central', ...
                       'FiniteDifferenceStepSize', 1e-6, ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance', 1e-10);

obj_fun = @(p) wls_error_m3(p, x0, u, DataFit, Probe_train);

% --- Multistart mit Latin Hypercube Sampling ---
N_lhs = 30;    %(Zum testen: 30) Anzahl LHS-Punkte fuer das billige Screening
K_opt = 2;      %(Zum testen: 2)  davon die besten -> teurer fmincon-Start
d     = numel(p0);

rng(42);                          % reproduzierbar
L = lhsdesign(N_lhs, d);          % N x d, platzfuellend in [0,1]

logLB = log10(pLB);  logUB = log10(pUB);
P = 10.^(logLB + L .* (logUB - logLB));   % N x d Startpunkte in [pLB,pUB]

% 1) Billiges Screening: Zielfunktion an allen LHS-Punkten (je 1 ODE-Lauf)
fprintf('LHS-Screening ueber %d Punkte ...\n', N_lhs);
Jscreen = inf(N_lhs,1);
for k = 1:N_lhs
    try, Jscreen(k) = obj_fun(P(k,:)); catch, end
end
[~, order] = sort(Jscreen); %Fehler sortieren - > Vielversprechendsten werden dann danach ausgewählt
Pstart = [p0; P(order(1:K_opt), :)];   % alter Startpunkt + beste K aus LHS

% 2) fmincon nur von den vielversprechendsten Startpunkten
p_opt = p0;  fval = inf;
for k = 1:size(Pstart,1)
    fprintf('--- fmincon Start %d/%d ---\n', k, size(Pstart,1));
    try
        [pk, Jk] = fmincon(obj_fun, Pstart(k,:), [], [], [], [], ...
                           pLB, pUB, [], options);
        if Jk < fval
            fval = Jk;  p_opt = pk;
        end
    catch Me
        fprintf('\nStart %d fehlgeschlagen:\n',k);
        fprintf('%s\n',Me.getReport);
    end
end

% 3. Ausgabe und Speichern
fprintf('\n--- Modell3 (ab Feed, t_feed=%.3f h): Fitting abgeschlossen ---\n', t_feed);
fprintf('WLS-Fehler (nur Punkte ab Feed): %.4f\n', fval);
% zum Vergleich: WLS desselben p_opt auf dem VOLLEN Datensatz
fval_full = wls_error_m3(p_opt, x0, u, DataFull, Probe_train);
fprintf('WLS-Fehler (voller Datensatz, zum Vergleich): %.4f\n', fval_full);
for i = 1:numel(p_opt)
    fprintf('%-10s = %.4f\n', namen{i}, p_opt(i));
end

scriptDir = pwd;
saveDir   = fullfile(scriptDir, '..', 'Daten');
save(fullfile(saveDir, 'p_opt_Modell3_woEtOH_fromfeed.mat'), 'p_opt', 't_feed');

%% Visualisierung -- Simulation ueber den VOLLEN Horizont, alle Punkte
t_start = u(1,1);
t_end   = max(DataFull.Biomasse.t(:)) + 1;
t_sim   = linspace(t_start, t_end, 300);

DOTstern = max(DataFull.O2.y);
options_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
X3 = sim_m3_sample(t_sim, x0, u, p_opt, DOTstern, Probe_train);

V    = X3(:,1);
cX   = X3(:,2)./V;   cGlc = X3(:,3)./V;   cAm = X3(:,4)./V; 4.
cPh  = X3(:,5)./V;   mB   = X3(:,6);      DOT = X3(:,7);

figure('Name','Modell3 - Fitting (ab Feed)','Position',[200 40 950 1000]);
plot_row(1, DataFull.Biomasse, t_sim, cX,   'Biomasse',  'c_X (g/L)',     t_feed, t_hi);
plot_row(2, DataFull.Glucose,  t_sim, cGlc, 'Glucose',   'c_{Glc} (g/L)', t_feed, t_hi);
plot_row(3, DataFull.Ammonium, t_sim, cAm,  'Ammonium',  'c_{Am} (g/L)',  t_feed, t_hi);
plot_row(4, DataFull.Phosphat, t_sim, cPh,  'Phosphat',  'c_{Ph} (g/L)',  t_feed, t_hi);
plot_row(5, DataFull.Base,     t_sim, mB,   'Base',      'm_B',           t_feed, t_hi);
plot_row(6, DataFull.O2,       t_sim, DOT,  'DOT',       'DOT (%)',       t_feed, t_hi);
xlabel('BatchAge (h)');


%% 5. Validierung auf ValData (RamScDef07) -- mit den auf TrainData gefitteten p_opt
u_val  = ValData.u;
x0_val = ValData.x0;
Probe_val = ValDataProbe;

% eigener Feed-Start des Validierungslaufs (gleiche Logik wie beim Training)
feedOn_val  = any(u_val(feedRows,:) > 0, 1);
idxFeed_val = find(feedOn_val, 1, 'first');
if isempty(idxFeed_val), t_feed_val = u_val(1,1); else, t_feed_val = u_val(1, idxFeed_val); end

% Validierungsfehler: nur Punkte ab Feed (vergleichbar zum Training)
ValFit   = window_data(ValData, t_feed_val, t_hi);
fval_val = wls_error_m3(p_opt, x0_val, u_val, ValFit, Probe_val);
fprintf('WLS-Fehler (Validierung, ab Feed t=%.3f h): %.4f\n', t_feed_val, fval_val);

% Simulation Validierung ueber den VOLLEN Horizont (ab x0 bei u(1,1))
t_end_val    = max(ValData.Biomasse.t(:)) + 1;
t_sim_val    = linspace(u_val(1,1), t_end_val, 300);
DOTstern_val = max(ValData.O2.y);
Xv = sim_m3_sample(t_sim_val, x0_val, u_val, p_opt, DOTstern_val, Probe_val);

Vv   = Xv(:,1);
cXv  = Xv(:,2)./Vv;   cGlcv = Xv(:,3)./Vv;   cAmv = Xv(:,4)./Vv;
cPhv = Xv(:,5)./Vv;   mBv   = Xv(:,6);       DOTv = Xv(:,7);

figure('Name','Modell3 - Validierung (RamScDef07)','Position',[220 40 950 1000]);
plot_row(1, ValData.Biomasse, t_sim_val, cXv,   'Biomasse (Val)',  'c_X (g/L)',     t_feed_val, t_hi);
plot_row(2, ValData.Glucose,  t_sim_val, cGlcv, 'Glucose (Val)',   'c_{Glc} (g/L)', t_feed_val, t_hi);
plot_row(3, ValData.Ammonium, t_sim_val, cAmv,  'Ammonium (Val)',  'c_{Am} (g/L)',  t_feed_val, t_hi);
plot_row(4, ValData.Phosphat, t_sim_val, cPhv,  'Phosphat (Val)',  'c_{Ph} (g/L)',  t_feed_val, t_hi);
plot_row(5, ValData.Base,     t_sim_val, mBv,   'Base (Val)',      'm_B',           t_feed_val, t_hi);
plot_row(6, ValData.O2,       t_sim_val, DOTv,  'DOT (Val)',       'DOT (%)',       t_feed_val, t_hi);
xlabel('BatchAge (h)');


% ======================= Hilfsfunktionen ==============================

function D = window_data(D, t_lo, t_hi)
% Behaelt aus jeder Messgroesse nur die Punkte mit t_lo <= t <= t_hi.
    flds = {'Biomasse','Glucose','Ammonium','Phosphat','Base','O2'};
    for i = 1:numel(flds)
        f = flds{i};
        m = D.(f);
        keep = (m.t >= t_lo) & (m.t <= t_hi);   % Orientierung bleibt erhalten
        m.t   = m.t(keep);
        m.y   = m.y(keep);
        m.var = m.var(keep);
        D.(f) = m;
    end
end

function report_window(Dfull, Dcut)
% Gibt aus, wie viele Punkte pro Messgroesse entfernt wurden.
    flds = {'Biomasse','Glucose','Ammonium','Phosphat','Base','O2'};
    for i = 1:numel(flds)
        f = flds{i};
        n0 = numel(Dfull.(f).t);
        n1 = numel(Dcut.(f).t);
        fprintf('  %-9s: %d -> %d  (%d entfernt)\n', f, n0, n1, n0-n1);
    end
end


%% WLS-Gütefunktional
function J = wls_error_m3(p, x0, u, Data, Probe)
% WLS für Modell3 ohne Ethanol.
% Zustände: x = [V; mX; mGlc; mAm; mPh; mB; DOT]
% Messgrössen: Biomasse->mX/V, Glucose->mGlc/V, Ammonium->mAm/V,
%               Phosphat->mPh/V, Base->mB, DOT->DOT.
    M = { Data.Biomasse, 'Biomasse',  2, true;   ...
          Data.Glucose,  'Glucose',   3, true;   ...
          Data.Ammonium, 'Ammonium',  4, true;   ...
          Data.Phosphat, 'Phosphat',  5, true;   ...
          Data.Base,     'Base',      6, false;  ...
          Data.O2,       'DOT',       7, false   };
    wmode = 'sum';
    wsig  = [1, 1, 1, 1, 1, 0.01]; %händische Gewichtungsmöglichkeit

    DOTstern = max(Data.O2.y);

    t_all = [];
    for i = 1:size(M,1)
        t_all = [t_all; M{i,1}.t(:)];
    end
    % Feed-Sprungzeiten zusätzlich aufnehmen
    tu = u(1,:).';
    tu = tu(tu >= min(t_all) & tu <= max(t_all));
    t_all = unique([t_all; tu]);
    try
    X = sim_m3_sample(t_all, x0, u, p, DOTstern, Probe);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:)))
        J = 1e8; return;
    end

    V = X(:,1);
    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1}; name = M{i,2}; idxState = M{i,3}; divByV = M{i,4};
        if isempty(mess.t), continue; end % Messreihe übersprungen, falls keine Daten vorhanden sind
        [tf, iT] = ismember(mess.t(:), t_all); % Prüfen, , ob alle Messzeitpunkte auch in den Simulationszeitpunkten vorkommen
        if any(~tf), warning('Messzeitpunkt für %s nicht gefunden.', name); J = 1e8; return; end
        y_sim = X(iT, idxState); % iT enthält Positionen der Messzeiten in t_all
        if divByV, y_sim = y_sim ./ V(iT); end
        % Formel für cost function aus Krämer 2017 Gl. 1
        r = (mess.y(:) - y_sim) ./ (max(mess.var(:), eps)); % max(), eps verhindert division durch 0
        switch wmode
            case 'sum',  contrib = sum(r.^2); % einfach aufsummeieren
            case 'mean', contrib = mean(r.^2); % noch 1/n teilen, damit unterschiedlich viele Messwerte ergebnis nicht verfälschen
        end
        J = J + wsig(i) * contrib;
    end
    if ~isfinite(J), J = 1e8; end
end



function plot_row(row, mess, t_sim, y_sim, name, ylab, t_feed, t_hi)
    subplot(7,1,row);
    errorbar(mess.t, mess.y, sqrt(mess.var), 'o', ...
             'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, y_sim, 'LineWidth', 2);
    xline(t_feed, '--k', 'Feed', 'LabelVerticalAlignment','middle');   % Feed-Start
    if nargin >= 8 && isfinite(t_hi)
        xline(t_hi, '--k', 't_{hi}', 'LabelVerticalAlignment','Top'); % Schwanz-Stop
    end
    title(['Modell3 - ' name]); ylabel(ylab);
    legend('Messung \pm \sigma','Simulation','Location','best'); grid on;
end



%% SAVE
%  Mit beiden bounds
%  --- Modell3 (ab Feed, t_feed=1.740 h): Fitting abgeschlossen ---
% WLS-Fehler (nur Punkte ab Feed): 1336.6468
% WLS-Fehler (voller Datensatz, zum Vergleich): 1299.1388
% mumax      = 0.1015
% KS         = 0.5400
% YXS        = 0.1207
% YAmX       = 0.3254
% YPhX       = 0.0376
% YB_Am      = 0.1036
% KLa        = 283.9009
% YXO        = 0.1657
% --- Modell3 (ab Feed, t_feed=1.740 h): Fitting abgeschlossen ---
% WLS-Fehler (nur Punkte ab Feed): 2281.2560
% WLS-Fehler (voller Datensatz, zum Vergleich): 2352.7695
% mumax      = 0.1456
% KS         = 0.8335
% YXS        = 0.3997
% YAmX       = 0.0963
% YPhX       = 0.0041
% YB_Am      = 0.1000
% KLa        = 292.5176
% YXO        = 0.9660
% --- Modell3 (ab Feed, t_feed=1.740 h): Fitting abgeschlossen ---
% WLS-Fehler (nur Punkte ab Feed): 1092.2166
% WLS-Fehler (voller Datensatz, zum Vergleich): 1049.4322 (110 / 4)
% mumax      = 0.1273
% KS         = 1.7301
% YXS        = 0.1660
% YAmX       = 0.2316
% YPhX       = 0.0601
% YB_Am      = 0.1000
% KLa        = 597.6215
% YXO        = 0.1004
% --- Modell3 (ab Feed, t_feed=1.740 h): Fitting abgeschlossen ---
% WLS-Fehler (nur Punkte ab Feed): 1028.9457
% WLS-Fehler (voller Datensatz, zum Vergleich): 898.4715??   (150 / 5)
% mumax      = 0.1617
% KS         = 4.9985
% YXS        = 0.1940
% YAmX       = 0.1835
% YPhX       = 0.0316
% YB_Am      = 0.1000
% KLa        = 281.7907
% YXO        = 0.2343
% --- Modell3 (ab Feed, t_feed=0.000 h): Fitting abgeschlossen ---
% WLS-Fehler (nur Punkte ab Feed): 1049.5533
% WLS-Fehler (voller Datensatz, zum Vergleich): 1049.5533       150 / 5 mit
% anfang und ende
% mumax      = 0.1479
% KS         = 2.7507
% YXS        = 0.2492
% YAmX       = 0.1255
% YPhX       = 0.1642
% YB_Am      = 1.1526
% KLa        = 32.1898
% YXO        = 3.0911
% --- Modell3 (ab Feed, t_feed=1.740 h): Fitting abgeschlossen ---
% WLS-Fehler (nur Punkte ab Feed): 1625.8064
% WLS-Fehler (voller Datensatz, zum Vergleich): 1499.6276       200 / 6
% mumax      = 0.0978
% KS         = 0.0920
% YXS        = 0.2067
% YAmX       = 0.1816
% YPhX       = 0.1188
% YB_Am      = 1.3420
% KLa        = 151.5621
% YXO        = 0.5454