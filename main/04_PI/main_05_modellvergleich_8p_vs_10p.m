%% main_05_modellvergleich_8p_vs_10p.m
%  Rechtfertigt die beiden Zusatzparameter KAm und KPh.
%
%  Verglichen werden auf IDENTISCHER Datenbasis und mit IDENTISCHEM
%  Guetefunktional:
%    Modell 3 alt (8p):  rX = mumax * cGlc/(KS+cGlc)
%    Modell 3 neu (10p): rX = mumax * cGlc/(KS+cGlc)
%                             * cAm/(KAm+cAm) * cPh/(KPh+cPh)
%
%  Damit der Vergleich kontrolliert ist, haben BEIDE Modelle dieselben
%  SECHS freien Parameter: mumax, YXS, YAmX, YPhX, YB_Am, YXO_eff.
%  KS und KLa sind in beiden fixiert (KS aus dem Batchfit Modell 1,
%  KLa als Reparametrisierung, s. main_04). Im 10p-Modell kommen KAm und
%  KPh hinzu, ebenfalls fixiert -- sie sind also KEINE zusaetzlichen
%  Freiheitsgrade fuer den Fit, sondern nur eine Strukturaenderung.
%  Ein besseres J kann deshalb nicht durch mehr Parameter erklaert werden.
%
%  VORAUSSETZUNG: in Modell3_woEtOH.m muss die Zeile
%      x = max(x, 0);
%  auskommentiert sein (wie im 10p-Modell). Sonst unterscheiden sich die
%  Modelle nicht nur in der Kinetik, sondern auch im numerischen Schutz.
% =========================================================================
clear; clc; close all;

projectRoot = pwd;
load(fullfile(projectRoot,'Daten/Daten_Processed/Processed_FedBatch_Modell3_MultiExp.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
addpath(fullfile(projectRoot,'..','utils'),'-begin');
rehash;

nTrain = numel(TrainSet);
nVal   = numel(ValSet);
DOTstern_train = arrayfun(@(D) max(D.O2.y), TrainSet);
DOTstern_val   = arrayfun(@(D) max(D.O2.y), ValSet);

wsig = [1 1 1 1 1 1];      % wie in main_04

fprintf('Training   : %s\n', strjoin({TrainSet.name}, ', '));
fprintf('Validierung: %s\n', strjoin({ValSet.name}, ', '));

%% 1. Parameter -----------------------------------------------------------
% Reihenfolge 8p : [mumax KS YXS YAmX YPhX YB_Am KLa YXO]
% Reihenfolge 10p: [... dieselben ... , KAm, KPh]
namen8 = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO'};

p0_8  = [ 0.30,  4.29,  0.30,  0.05,  0.02,  0.02,  1.0,  200 ];
pLB_8 = [ 0.01,  4.29,  0.05,  1e-3,  1e-4,  1e-4,  1.0,  0.1 ];
pUB_8 = [ 1.00,  4.29,  1.00,  1.00,  1.00,  1.0,   1.0,  1e5 ];

iFix8 = [2 7];                 % KS, KLa  -> pLB == pUB haelt sie fest
isFree8 = true(1,8);  isFree8(iFix8) = false;

% Ergebnis des 10p-Fits als Startpunkt (die ersten 8 Eintraege passen 1:1)
p_opt_10 = [0.3193 4.29 0.1309 0.0961 0.0926 0.0218 1.0 209.5958 0.01 0.01];
pRef_8   = p_opt_10(1:8);
pRef_8   = min(max(pRef_8, pLB_8), pUB_8);

options = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                       'MaxFunctionEvaluations', 8000, ...
                       'FiniteDifferenceType','central', ...
                       'FiniteDifferenceStepSize', 1e-6, ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance', 1e-10);

MODELL8  = @Modell3_woEtOH;        % altes Modell, einfache Monod-Kinetik
obj8 = @(p) wls_multi(p, TrainSet, DOTstern_train, wsig, MODELL8);

%% 2. Fit des alten Modells ----------------------------------------------
J_ref8 = obj8(pRef_8);
fprintf('\nStartpunkt (aus 10p-Optimum): J = %.4f\n', J_ref8);

N_lhs = 15;  K_opt = 1;
rng(42);
L = lhsdesign(N_lhs, 8);
P = 10.^(log10(pLB_8) + L .* (log10(pUB_8) - log10(pLB_8)));
P(:,iFix8) = repmat(pLB_8(iFix8), N_lhs, 1);

fprintf('LHS-Screening ueber %d Punkte ...\n', N_lhs);
Js = inf(N_lhs,1);
for k = 1:N_lhs
    try, Js(k) = obj8(P(k,:)); catch, end
end
[~, ord] = sort(Js);
Pstart = [pRef_8; p0_8; P(ord(1:min(K_opt,sum(isfinite(Js)))), :)];

p_opt_8 = pRef_8;  fval_8 = J_ref8;
for k = 1:size(Pstart,1)
    fprintf('\n--- fmincon Start %d/%d ---\n', k, size(Pstart,1));
    try
        [pk, Jk] = fmincon(obj8, Pstart(k,:), [], [], [], [], ...
                           pLB_8, pUB_8, [], options);
        if Jk < fval_8, fval_8 = Jk;  p_opt_8 = pk; end
    catch ME
        fprintf('  fehlgeschlagen: %s\n', ME.message);
    end
end

%% 3. Ergebnis des alten Modells -----------------------------------------
fprintf('\n=====================================================\n');
fprintf('  Modell 3 ALT (8p, einfache Monod-Kinetik)\n');
fprintf('=====================================================\n');
fprintf('WLS gesamt (Training): %.4f\n\n', fval_8);
for i = 1:8
    if isFree8(i)
        fprintf('%-8s = %10.4f\n', namen8{i}, p_opt_8(i));
    else
        fprintf('%-8s = %10.4f   (fixiert)\n', namen8{i}, p_opt_8(i));
    end
end

%% 4. Direkter Vergleich --------------------------------------------------
MODELL10 = @Modell3_woEtOH_10p;

fprintf('\n=====================================================\n');
fprintf('  VERGLEICH  (beide: 6 freie Parameter, gleiche Daten)\n');
fprintf('=====================================================\n');

[Jtr8,  ch8]  = eval_set(p_opt_8,  TrainSet, DOTstern_train, wsig, MODELL8);
[Jtr10, ch10] = eval_set(p_opt_10, TrainSet, DOTstern_train, wsig, MODELL10);
[Jv8,   cv8]  = eval_set(p_opt_8,  ValSet,   DOTstern_val,   wsig, MODELL8);
[Jv10,  cv10] = eval_set(p_opt_10, ValSet,   DOTstern_val,   wsig, MODELL10);

fprintf('\n%-14s %12s %12s %10s\n', '', 'alt (8p)', 'neu (10p)', 'Aenderung');
fprintf('%-14s %12.1f %12.1f %9.1f %%\n', 'J Training', sum(Jtr8), sum(Jtr10), ...
        100*(sum(Jtr10)-sum(Jtr8))/sum(Jtr8));
fprintf('%-14s %12.1f %12.1f %9.1f %%\n', 'J Validierung', sum(Jv8), sum(Jv10), ...
        100*(sum(Jv10)-sum(Jv8))/sum(Jv8));

fprintf('\n--- Beitrag pro Experiment (Training) ---\n');
fprintf('%-14s %12s %12s\n', 'Experiment', 'alt (8p)', 'neu (10p)');
for k = 1:nTrain
    fprintf('%-14s %12.1f %12.1f\n', TrainSet(k).name, Jtr8(k), Jtr10(k));
end

kn = {'Biomasse','Glucose','Ammonium','Phosphat','Base','DOT'};
fprintf('\n--- Beitrag pro Kanal (Summe Training) ---\n');
fprintf('%-14s %12s %12s\n', 'Kanal', 'alt (8p)', 'neu (10p)');
for i = 1:6
    fprintf('%-14s %12.1f %12.1f\n', kn{i}, sum(ch8(:,i)), sum(ch10(:,i)));
end

fprintf('\n--- Beitrag pro Kanal (Validierung) ---\n');
fprintf('%-14s %12s %12s\n', 'Kanal', 'alt (8p)', 'neu (10p)');
for i = 1:6
    fprintf('%-14s %12.1f %12.1f\n', kn{i}, sum(cv8(:,i)), sum(cv10(:,i)));
end

%% 5. Das eigentliche Argument: negative Zustaende ------------------------
% Ohne die Monod-Faktoren laeuft rX auch dann weiter, wenn Ammonium oder
% Phosphat aufgebraucht sind -> die Massenbilanz zieht mAm/mPh unter null.
% Das ist physikalisch unmoeglich und nicht bloss ein Fit-Problem.
% Deshalb hier OHNE die NonNegative-Option integrieren: nur so wird
% sichtbar, was das Modell selbst tut, statt was der Integrator abfaengt.
fprintf('\n=====================================================\n');
fprintf('  NEGATIVE ZUSTAENDE (Integration ohne NonNegative)\n');
fprintf('=====================================================\n');
fprintf('%-14s %10s %10s %10s %10s\n', 'Experiment', ...
        'min mAm 8p', 'min mPh 8p', 'min mAm 10p', 'min mPh 10p');

Alle = [TrainSet, ValSet];
DSa  = [DOTstern_train, DOTstern_val];
for k = 1:numel(Alle)
    D = Alle(k);
    t = sort(unique([D.Biomasse.t; D.Ammonium.t; D.Phosphat.t]));
    a = min_states(p_opt_8,  D, DSa(k), MODELL8,  t);
    b = min_states(p_opt_10, D, DSa(k), MODELL10, t);
    fprintf('%-14s %10.3f %10.3f %10.3f %10.3f\n', D.name, a(1), a(2), b(1), b(2));
end
fprintf(['\nNegative Werte bei 8p belegen: die Zusatzterme sind kein\n', ...
         'Fit-Trick, sondern stellen die physikalische Konsistenz her.\n']);


%% ======================================================================
%  Lokale Funktionen
%  ======================================================================

function [Jv, contrib] = eval_set(p, Set, DOTstern, wsig, MODELFUN)
% Wertet jedes Experiment einzeln aus -> J je Experiment + Kanalbeitraege.
    n = numel(Set);
    Jv = zeros(n,1);  contrib = zeros(n,6);
    for k = 1:n
        [Jv(k), contrib(k,:)] = wls_one(p, Set(k).x0, Set(k).u, Set(k), ...
                                        Set(k).Probe, DOTstern(k), wsig, MODELFUN);
    end
end


function J = wls_multi(p, Set, DOTstern, wsig, MODELFUN)
    J = 0;
    for k = 1:numel(Set)
        Jk = wls_one(p, Set(k).x0, Set(k).u, Set(k), Set(k).Probe, ...
                     DOTstern(k), wsig, MODELFUN);
        J = J + Jk;
        if ~isfinite(J), J = 1e8; return; end
    end
end


function [J, contrib, nch] = wls_one(p, x0, u, Data, Probe, DOTstern, wsig, MODELFUN)
% Identisch zu wls_error_m3 in main_04 -- nur das Modell ist austauschbar.
% Trapez-Zeitgewichte pro Kanal (Summe 1), Base als Inkrement.
    M = { Data.Biomasse, 2, true;   Data.Glucose,  3, true; ...
          Data.Ammonium, 4, true;   Data.Phosphat, 5, true; ...
          Data.Base,     6, false;  Data.O2,       7, false };
    nm = {'Biomasse','Glucose','Ammonium','Phosphat','Base','DOT'};

    contrib = zeros(1,6);  nch = 0;

    t_all = [];
    for i = 1:6
        if wsig(i) == 0, continue; end
        t_all = [t_all; M{i,1}.t(:)]; %#ok<AGROW>
    end
    if isempty(t_all), J = 1e8; return; end
    tu = u(1,:).';
    t_all = unique([t_all; tu(tu >= min(t_all) & tu <= max(t_all))]);

    ws = warning('off','MATLAB:ode15s:IntegrationTolNotMet');
    try
        X = sim_generic(t_all, x0, u, p, DOTstern, Probe, MODELFUN, true);
    catch
        warning(ws);  J = 1e8;  return;
    end
    warning(ws);
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:))), J = 1e8; return; end

    V = X(:,1);  J = 0;
    for i = 1:6
        if wsig(i) == 0, continue; end
        mess = M{i,1};  idxState = M{i,2};  divByV = M{i,3};
        if isempty(mess.t), continue; end

        [tf, iT] = ismember(mess.t(:), t_all);
        if any(~tf), J = 1e8; return; end
        y_sim = X(iT, idxState);
        if divByV, y_sim = y_sim ./ V(iT); end

        if strcmp(nm{i},'Base')
            keep = subsample_idx(mess.t, 2.0);
            if nnz(keep) < 2, continue; end
            tB = mess.t(keep);  yB = mess.y(keep);
            vB = max(mess.var(keep), eps);  ysB = y_sim(keep);
            r = (diff(yB) - diff(ysB)) ./ sqrt(vB(2:end) + vB(1:end-1));
            w = time_weights(tB(2:end));
        else
            r = (mess.y(:) - y_sim) ./ sqrt(max(mess.var(:), eps));
            w = time_weights(mess.t);
        end

        contrib(i) = wsig(i) * sum(w .* r.^2);
        J          = J + contrib(i);
        nch        = nch + 1;
    end
    if ~isfinite(J) || J < 0, J = 1e8; end
end


function X = sim_generic(t, x0, u, p, DOTstern, Probe, MODELFUN, nonNeg)
% Wie sim_m3_sample_10p, aber mit austauschbarem Modell und schaltbarem
% NonNegative. Segmentgrenzen nur an echten Unstetigkeiten (Feed, Probe).
    t  = t(:);
    tp = Probe.BatchAge(:);  vp = Probe.Volumen(:);
    tf = u(1,:).';           t0 = u(1,1);

    I  = tp >= t0 & tp <= max(t);   tp = tp(I);  vp = vp(I);
    tf = tf(tf >= t0 & tf <= max(t));
    tb = unique([t0; tf; tp; max(t)]);

    x = x0(:);  X = zeros(numel(t), numel(x0));

    if nonNeg
        opt = odeset('RelTol',1e-8,'AbsTol',1e-10,'NonNegative',1:7,'MaxStep',1.0);
    else
        opt = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1.0);
    end

    X(t == t0,:) = repmat(x.', nnz(t == t0), 1);
    for j = find(tp == t0).', x = probe_m3(x, vp(j)); end

    for k = 2:numel(tb)
        ta = tb(k-1);  te = tb(k);
        sel  = t > ta & t <= te;
        tout = unique([ta; t(sel); te]);

        [~, Xs] = ode15s(@(tt,xx) MODELFUN(tt,xx,u,p,DOTstern), tout, x, opt);
        if numel(tout) == 2, Xs = Xs([1 end],:); end

        for i = find(sel).'
            X(i,:) = Xs(find(tout == t(i), 1), :);
        end
        x = Xs(end,:).';
        for j = find(tp == te).', x = probe_m3(x, vp(j)); end
    end
end


function mn = min_states(p, D, DOTstern, MODELFUN, t)
% Kleinste erreichte Ammonium-/Phosphat-Masse OHNE NonNegative.
    try
        X = sim_generic(t, D.x0, D.u, p, DOTstern, D.Probe, MODELFUN, false);
        mn = [min(X(:,4)), min(X(:,5))];
    catch
        mn = [NaN NaN];       % Integration abgebrochen -> ebenfalls ein Befund
    end
end


function w = time_weights(t)
% Trapez-Zeitgewichte, normiert auf Summe 1.
    t = t(:);  n = numel(t);
    if n == 0, w = [];  return; end
    if n == 1, w = 1;   return; end
    a = zeros(n,1);
    a(1) = t(2) - t(1);  a(end) = t(end) - t(end-1);
    a(2:end-1) = (t(3:end) - t(1:end-2)) / 2;
    T = t(end) - t(1);
    if T <= 0, w = ones(n,1)/n; return; end
    w = max(a, 0) / T;
end


function keep = subsample_idx(t, dt_min)
    keep = false(numel(t),1);  last = -inf;
    for i = 1:numel(t)
        if t(i) - last >= dt_min, keep(i) = true;  last = t(i); end
    end
end
