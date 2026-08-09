% main_06_ovp_modell3.m
%
% Optimale Versuchsplanung (OVP) fuer Modell 3.
%
% Frage: wie muss der Glucose-Feed eines ZUSAETZLICHEN Experiments
% verlaufen, damit die Parameter danach moeglichst genau bekannt sind?
% Die Parameter selbst bleiben dabei fest -- optimiert wird nur u_Glc(t).
%
% Vorgehen:
%   1. FM_alt aus den vier Trainingsexperimenten (aus main_05).
%   2. Fuer einen Kandidaten-Feedverlauf das neue Experiment simulieren
%      (Zustaende + Sensitivitaeten) und daraus FM_neu bilden. Die
%      Messvarianz kommt aus sigma^2 = a*y + b, ausgewertet auf der
%      SIMULIERTEN Groesse -- Messwerte gibt es ja noch keine.
%   3. A-Kriterium auf der NORMIERten Gesamt-FIM minimieren:
%        J = trace( inv( diag(p) * (FM_alt + FM_neu) * diag(p) ) )
%      also die Summe der relativen Parametervarianzen. Ohne Normierung
%      wuerde YXO (~246) gegenueber YB_Am (~0.02) alles dominieren.
%   4. Nebenbedingungen: Feedgrenzen, Reaktorvolumen, Glucose, DOT.
%
% Nur die sechs freien Parameter gehen ein -- ueber alle zehn waere die
% FIM singulaer.

clear; clc; close all;

projectRoot = pwd;
DATEN = fullfile(projectRoot,'Daten');        % zentrale Pfadwurzel
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
addpath(fullfile(projectRoot,'..','utils'),  '-begin');
rehash;

%% 1. Eingaenge laden -----------------------------------------------------
load(fullfile(DATEN,'Daten_Processed','Processed_FedBatch_Modell3_MultiExp.mat'));
S  = load(fullfile(DATEN,'p_opt','p_opt_Modell3_woEtOH_10p_multi.mat'));
F0 = load(fullfile(DATEN,'FIM','FIM_Modell3_multi.mat'));

p     = S.p_opt(:);
namen = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO','KAm','KPh'};
iFree = [1 3 4 5 6 8];
nf    = numel(iFree);
nx    = 7;   np = 10;

FM_alt = F0.FM;          % 10x10, aus den vier Trainingsexperimenten

% FM_alt und p muessen aus DEMSELBEN Fit stammen, sonst beschreibt die
% Ausgangslage einen anderen Punkt als die Optimierung.
if max(abs(F0.p_opt(:) - p)) > 1e-9
    error(['FIM_Modell3_multi.mat wurde mit anderen Parametern gerechnet. ', ...
           'Erst main_05_uncertainty_FIM_m3_multi.m neu laufen lassen.']);
end

% Varianzparameter a,b aus einem Rohdatensatz (gleiche Quelle wie im
% Preprocessing), Einheitenumrechnung wie dort.
Mess = load(fullfile(DATEN,'MessDaten_SPI1_Projekt','Mess_RamScDef10.mat')).Mess;
VP.ab = [Mess.Messdaten.Biomasse.VarParam.a Mess.Messdaten.Biomasse.VarParam.b;
         Mess.Messdaten.Glucose.VarParam.a  Mess.Messdaten.Glucose.VarParam.b;
         Mess.Messdaten.Ammonium.VarParam.a Mess.Messdaten.Ammonium.VarParam.b;
         Mess.Messdaten.Phosphat.VarParam.a Mess.Messdaten.Phosphat.VarParam.b;
         Mess.Messdaten.BASE.VarParam.a     Mess.Messdaten.BASE.VarParam.b;
         Mess.Messdaten.O2.VarParam.a       Mess.Messdaten.O2.VarParam.b];

fprintf('Freie Parameter: %s\n', strjoin(namen(iFree), ', '));

%% 2. Auslegung des geplanten Experiments --------------------------------
VP.T      = 15;        % Horizont [h]
VP.nu     = 10;        % stueckweise konstante Feed-Abschnitte
VP.dt_off = 1.5;       % Abtastung Offline (Biomasse, Glc, Am, Ph) [h]
VP.dt_on  = 0.1;       % Abtastung Online (Base, DOT) [h]
VP.Vprobe = 0.03;      % Probenvolumen je Offline-Probe [L]

VP.p        = p;
VP.DOTstern = max(TrainSet(2).O2.y);    % wie im Fit, nicht von Hand gesetzt
VP.x0       = TrainSet(2).x0;           % Startzustand wie RamScDef04
VP.x0(1)    = 10;                       % V0 = 10 L
VP.x0(3)    = 3 * VP.x0(1);             % cGlc0 ~ 3 g/L, also unter KS --
                                        % sonst verbringt der Lauf den
                                        % halben Horizont damit, die
                                        % Anfangscharge abzubauen

% Feste Feeds und Feedkonzentrationen
VP.uAm = 0.005;   VP.cAm_in = 30;
VP.uPh = 0.002;   VP.cPh_in = 24;
VP.cGlc_in = 450;

VP.tu    = linspace(0, VP.T, VP.nu+1);   % Stuetzstellen des Feedverlaufs
VP.t_off = (0:VP.dt_off:VP.T).';
VP.t_on  = (0:VP.dt_on :VP.T).';

% Probenahmen an den Offline-Zeitpunkten (ausser t = 0)
VP.Probe.BatchAge = VP.t_off(2:end);
VP.Probe.Volumen  = VP.Vprobe * ones(numel(VP.t_off)-1, 1);

%% 3. Nebenbedingungen ----------------------------------------------------
VP.uGlc_max = 0.20;    % L/h
VP.Vmax     = 14.0;    % L, Reaktorgrenze
VP.cGlc_max = 35;      % g/L, aus vergangenen Versuchen
VP.DOTmin   = 20;      % %, keine O2-Limitierung erzwingen

%% 4. Ausgangslage: was wissen wir ohne neues Experiment? ----------------
[relStd_alt, CN_alt, CV_alt] = analyse_fim(FM_alt, p, iFree);

fprintf('\n--- VOR OVP (nur Altversuche) ---\n');
fprintf('%-9s %12s %10s\n','Param','Wert','rel. [%]');
for i = 1:nf
    fprintf('%-9s %12.4f %10.1f\n', namen{iFree(i)}, p(iFree(i)), 100*relStd_alt(i));
end
fprintf('A-Kriterium = %.4e | Kondition = %.3e\n', trace(CV_alt), CN_alt);

%% 5. Feedverlauf optimieren ---------------------------------------------
u0  = 0.05 * ones(1, VP.nu);
uLB = zeros(1, VP.nu);
uUB = VP.uGlc_max * ones(1, VP.nu);

obj = @(uq) ovp_criterion(uq, VP, FM_alt, iFree, nx, np);
con = @(uq) ovp_constraints(uq, VP, nx);

opts = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 3000, ...
                    'FiniteDifferenceStepSize', 1e-4);

fprintf('\nStarte OVP (%d Abschnitte, Horizont %.0f h) ...\n', VP.nu, VP.T);
J0 = obj(u0);
fprintf('Startwert (konstanter Feed): A = %.4e\n', J0);

[u_opt, J_opt] = fmincon(obj, u0, [], [], [], [], uLB, uUB, con, opts);

%% 6. Ergebnis ------------------------------------------------------------
[FM_neu, X, t_sim] = fim_design(u_opt, VP, nx, np);
FM_ges = FM_alt + FM_neu;
[relStd_neu, CN_neu, ~] = analyse_fim(FM_ges, p, iFree);

% Drei Vergleichsstufen: ohne neues Experiment, mit konstantem Feed, mit
% optimiertem Feed. Der Sprung von "ohne" auf "konstant" ist der Nutzen
% des Experiments an sich, der Rest ist der Nutzen der Optimierung.
fprintf('\n--- A-Kriterium ---\n');
fprintf('nur Altversuche      : %.4e\n', trace(CV_alt));
fprintf('+ konstanter Feed    : %.4e\n', J0);
fprintf('+ optimierter Feed   : %.4e   (%.1f %% besser als konstant)\n', ...
        J_opt, 100*(J0-J_opt)/J0);

fprintf('\n--- Parameterunsicherheiten ---\n');
fprintf('%-9s %12s %12s %12s %10s\n', ...
        'Param','Wert','rel. alt [%]','rel. neu [%]','Faktor');
for i = 1:nf
    fprintf('%-9s %12.4f %12.1f %12.1f %10.2f\n', namen{iFree(i)}, p(iFree(i)), ...
            100*relStd_alt(i), 100*relStd_neu(i), relStd_alt(i)/relStd_neu(i));
end
fprintf('Kondition: %.3e -> %.3e\n', CN_alt, CN_neu);

save(fullfile(DATEN,'FIM','OVP_Modell3.mat'), ...
     'u_opt','J_opt','J0','FM_neu','FM_ges','relStd_alt','relStd_neu', ...
     'CN_alt','CN_neu','VP','iFree');

%% 7. Abbildungen ---------------------------------------------------------
V = X(:,1);
figure('Name','OVP Modell 3','Position',[150 60 900 900]);

subplot(4,1,1);
stairs(VP.tu, [u_opt u_opt(end)], 'LineWidth', 2); hold on;
stairs(VP.tu, [u0 u0(end)], '--', 'LineWidth', 1.2);
ylabel('u_{Glc} [L/h]'); title('Optimaler Glucose-Feed');
legend('optimiert','Startwert','Location','best'); grid on;

subplot(4,1,2);
plot(t_sim, X(:,3)./V, 'LineWidth', 2); hold on;
yline(p(2), '--r', 'K_S');
ylabel('c_{Glc} [g/L]'); title('Glucose'); grid on;

subplot(4,1,3);
plot(t_sim, X(:,2)./V, 'LineWidth', 2); hold on;
plot(t_sim, X(:,4)./V, 'LineWidth', 1.5);
plot(t_sim, X(:,5)./V, 'LineWidth', 1.5);
ylabel('c [g/L]'); legend('c_X','c_{Am}','c_{Ph}','Location','best'); grid on;

subplot(4,1,4);
yyaxis left;  plot(t_sim, V, 'LineWidth', 2);        ylabel('V [L]');
yyaxis right; plot(t_sim, X(:,7), 'LineWidth', 2);   ylabel('DOT [%]');
xlabel('Zeit [h]'); grid on;

figure('Name','Parameterunsicherheit vor/nach OVP');
bar(100*[relStd_alt(:) relStd_neu(:)]);
set(gca,'XTickLabel', namen(iFree));
ylabel('relative Standardabweichung [%]');
legend('vor OVP','nach OVP','Location','best');
title(sprintf('A-Kriterium %.2e -> %.2e', J0, J_opt)); grid on;


%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function J = ovp_criterion(uq, VP, FM_alt, iFree, nx, np)
% A-Kriterium: Summe der relativen Parametervarianzen -> minimieren.
% Alternativen waeren D (det) oder modifiziertes E (max/min Eigenwert).
% Das A-Kriterium verbessert bevorzugt die ohnehin gut bestimmten
% Richtungen und laesst die schlechteste weitgehend unberuehrt.
    try
        FM_neu = fim_design(uq, VP, nx, np);
    catch
        J = 1e12; return;                 % Integration gescheitert
    end
    Ff = FM_alt(iFree,iFree) + FM_neu(iFree,iFree);
    pf = VP.p(iFree);
    Fn = diag(pf) * Ff * diag(pf);        % dimensionslos

    if rcond(Fn) < 1e-14 || any(~isfinite(Fn(:)))
        J = 1e12; return;
    end
    J = trace(Fn \ eye(numel(iFree)));
    if ~isfinite(J) || J <= 0, J = 1e12; end
end


function [c, ceq] = ovp_constraints(uq, VP, nx)
% Zustandsbeschraenkungen. Hier reicht die Simulation ohne
% Sensitivitaeten -- das ist deutlich billiger als die volle FIM.
    ceq = [];
    try
        X = sim_states(uq, VP);
    catch
        c = 1e3; return;
    end
    V = X(:,1);
    c = [ max(V) - VP.Vmax;                   % V   <= Vmax
          max(X(:,3)./V) - VP.cGlc_max;       % cGlc <= cGlc_max
          VP.DOTmin - min(X(:,7)) ];          % DOT >= DOTmin
end


function [FM, X, t_sim] = fim_design(uq, VP, nx, np)
% Simuliert das geplante Experiment mit Sensitivitaeten und bildet die FIM.
    u     = build_u(uq, VP);
    t_sim = unique([VP.t_off; VP.t_on; VP.tu(:)]);

    X_ext = sim_core(t_sim, [VP.x0(:); zeros(nx*np,1)], u, VP.p, ...
                     VP.DOTstern, VP.Probe, @Modell3_woEtOH_XP, nx, np);
    X = X_ext(:, 1:nx);

    % {Index in der Messgleichung, Zeitraster, Zeile in VP.ab}
    kan = { 1, VP.t_off, 1;    % Biomasse
            2, VP.t_off, 2;    % Glucose
            3, VP.t_off, 3;    % Ammonium
            4, VP.t_off, 4;    % Phosphat
            5, VP.t_on,  5;    % Base
            6, VP.t_on,  6 };  % DOT

    FM = zeros(np, np);
    for i = 1:6
        iOut = kan{i,1};
        [tf, idx] = ismember(kan{i,2}, t_sim);
        if any(~tf), error('Zeitpunkt nicht im Simulationsraster.'); end

        n  = numel(idx);
        Sy = zeros(n, np);
        yv = zeros(n, 1);
        for k = 1:n
            xk   = X_ext(idx(k), 1:nx).';
            XP_k = reshape(X_ext(idx(k), nx+1:end), nx, np);
            dhdx    = Modell3_dmgldx(xk);              % ny x nx
            Sy(k,:) = dhdx(iOut,:) * XP_k;             % Kettenregel
            yk      = Modell3_mgl(xk);
            yv(k)   = yk(iOut);
        end

        % sigma^2 = a*y + b auf der SIMULIERTEN Groesse
        a = VP.ab(kan{i,3},1);   b = VP.ab(kan{i,3},2);
        switch iOut
            case 5, v = a.*yv/1000 + b/1e6;     % Base: mL -> L
            case 6, v = 100*a.*yv + 1e4*b;      % O2: Anteil -> %
            otherwise, v = a.*yv + b;
        end
        v = max(v, eps);

        if iOut == 5
            % Base als INKREMENT, exakt wie in FM_alt (main_05). Sonst
            % werden zwei Matrizen mit unterschiedlicher Konvention addiert.
            keep = subsample_idx(kan{i,2}, 2.0);
            if nnz(keep) < 3, continue; end
            dS = diff(Sy(keep,:), 1, 1);
            vk = v(keep);
            vd = vk(2:end) + vk(1:end-1);
            m  = size(dS,1);
            for k = 1:m
                FM = FM + (dS(k,:).' * dS(k,:)) / vd(k) / m;
            end
        else
            % Mittelung pro Kanal -- gleiche Konvention wie im Guetefunktional
            for k = 1:n
                FM = FM + (Sy(k,:).' * Sy(k,:)) / v(k) / n;
            end
        end
    end
end




function X = sim_states(uq, VP)
% Nur Zustaende, ohne Sensitivitaeten (fuer die Nebenbedingungen).
    u     = build_u(uq, VP);
    t_sim = unique([VP.t_off; VP.t_on; VP.tu(:)]);
    X = sim_core(t_sim, VP.x0(:), u, VP.p, VP.DOTstern, VP.Probe, ...
                 @Modell3_woEtOH_10p, 7, 0);
end


function u = build_u(uq, VP)
% Baut die u-Matrix aus dem stueckweise konstanten Glucose-Feed.
    u = zeros(10, numel(VP.tu));
    u(1,:) = VP.tu;
    u(2,:) = VP.uAm;   u(3,:) = VP.cAm_in;
    u(4,:) = VP.uPh;   u(5,:) = VP.cPh_in;
    u(6,:) = [uq(:).' uq(end)];        % letzter Wert wird gehalten
    u(7,:) = VP.cGlc_in;
end


function X = sim_core(t, x0, u, p, DOTstern, Probe, rhs, nx, np)
% Segmentweise Integration mit Probenahme-Spruengen. Segmentgrenzen nur an
% echten Unstetigkeiten (Feedwechsel, Probenahme) -- sonst haengt das
% Ergebnis vom Auswerteraster ab.
    t  = t(:);
    tp = Probe.BatchAge(:);  vp = Probe.Volumen(:);
    tf = u(1,:).';           t0 = u(1,1);

    I  = tp >= t0 & tp <= max(t);  tp = tp(I);  vp = vp(I);
    tf = tf(tf >= t0 & tf <= max(t));
    tb = unique([t0; tf; tp; max(t)]);

    x = x0(:);
    X = zeros(numel(t), numel(x0));
    opt = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);

    X(t == t0,:) = repmat(x.', nnz(t == t0), 1);
    for j = find(tp == t0).', x = probe_any(x, vp(j), nx, np); end

    for k = 2:numel(tb)
        ta = tb(k-1);  te = tb(k);
        sel  = t > ta & t <= te;
        tout = unique([ta; t(sel); te]);

        [~, Xs] = ode15s(@(tt,xx) rhs(tt,xx,u,p,DOTstern), tout, x, opt);
        if numel(tout) == 2, Xs = Xs([1 end],:); end
        if size(Xs,1) ~= numel(tout), error('Integration unvollstaendig.'); end

        for i = find(sel).'
            X(i,:) = Xs(find(tout == t(i), 1), :);
        end
        x = Xs(end,:).';
        for j = find(tp == te).', x = probe_any(x, vp(j), nx, np); end
    end
end


function xe = probe_any(xe, Vp, nx, np)
% Probenahme auf den Zustand und -- falls vorhanden -- auf die
% Sensitivitaeten anwenden: x+ = g(x), also XP+ = (dg/dx) * XP.
    x = xe(1:nx);
    if np == 0
        xe = probe_m3(x, Vp);
        return
    end
    XP = reshape(xe(nx+1:end), nx, np);

    Jg = zeros(nx);                       % dg/dx numerisch
    for j = 1:nx
        h  = 1e-6 * max(abs(x(j)), 1);
        xp = x; xp(j) = xp(j) + h;
        xm = x; xm(j) = xm(j) - h;
        Jg(:,j) = (probe_m3(xp,Vp) - probe_m3(xm,Vp)) / (2*h);
    end

    xe = [reshape(probe_m3(x,Vp), [], 1); reshape(Jg*XP, [], 1)];
end


function [relStd, CN, CV] = analyse_fim(FM, p, iFree)
% Relative Standardabweichungen und Konditionszahl der NORMIERten FIM.
% Durch die Normierung sind die Diagonalelemente von CV bereits relative
% Varianzen.
    pf = p(iFree);
    Fn = diag(pf) * FM(iFree,iFree) * diag(pf);

    if rcond(Fn) < 1e-14
        CV = pinv(Fn);
    else
        CV = Fn \ eye(numel(iFree));
    end
    relStd = sqrt(abs(diag(CV)));
    en = sort(eig(Fn),'descend');
    CN = en(1)/en(end);
end

    function keep = subsample_idx(t, dt_min)
        % Waehlt Punkte mit mindestens dt_min Abstand (erster Punkt immer dabei).
        keep = false(numel(t),1);  last = -inf;
        for i = 1:numel(t)
            if t(i) - last >= dt_min, keep(i) = true;  last = t(i); end
        end
    end