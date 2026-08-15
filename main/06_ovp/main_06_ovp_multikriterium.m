% main_06_ovp_multikriterium.m
%
% Optimale Versuchsplanung (OVP) fuer Modell 3.
%
% Entscheidungsvariablen: Glucose-, Ammonium- und Phosphatfeed, normiert
% auf [0,1]. Am/Ph laufen auf einem groeberen Raster als Glucose, weil ihr
% informativer Gehalt ein Umschalten ist und keine feine Form.
%
% Kriterien: A (Summe der Varianzen) und D (Volumen des Ellipsoids).
% D ist invariant gegenueber linearer Reparametrisierung -- die
% diag(p)-Normierung verschiebt D nur um eine Konstante und aendert das
% Optimum nicht. Fuer A ist die Normierung dagegen zwingend, sonst
% dominiert YXO (~0.5) gegenueber YB_Am (~0.03).
%
% Die FIM muss denselben Schaetzer beschreiben wie der Fit:
%   - Gewichte wsig aus p_opt-Datei
%   - Mittelung pro Kanal (Faktor 1/n)
%   - Base als Absolutwert (Kraemer & King 2017, Gl. 22)
%   - keine O2-Umrechnung mehr (a,b sind bereits in %)

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fullfile(projectRoot,'..','..');
DATEN = fullfile(projectRoot,'Daten');
addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');
rehash;

%% 1. Eingaenge laden -----------------------------------------------------
load(fullfile(DATEN,'Daten_Processed','Processed_FedBatch_Modell3_MultiExp.mat'));
S  = load(fullfile(DATEN,'p_opt','p_opt_Modell3_woEtOH_10p_multi.mat'));
F0 = load(fullfile(DATEN,'FIM','FIM_Modell3_multi.mat'));

p     = S.p_opt(:);
namen = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO','KAm','KPh'};
iFree = [1 2 3 4 5 6 8];
nf    = numel(iFree);
FM_alt = F0.FM;

if max(abs(F0.p_opt(:) - p)) > 1e-9
    error('FIM-Datei passt nicht zu p_opt -- erst main_05 neu laufen lassen.');
end

% Varianzparameter a,b (gleiche Quelle wie im Preprocessing)
Mess = load(fullfile(DATEN,'MessDaten_SPI1_Projekt','Mess_RamScDef10.mat')).Mess;
VP.ab = [Mess.Messdaten.Biomasse.VarParam.a Mess.Messdaten.Biomasse.VarParam.b;
         Mess.Messdaten.Glucose.VarParam.a  Mess.Messdaten.Glucose.VarParam.b;
         Mess.Messdaten.Ammonium.VarParam.a Mess.Messdaten.Ammonium.VarParam.b;
         Mess.Messdaten.Phosphat.VarParam.a Mess.Messdaten.Phosphat.VarParam.b;
         Mess.Messdaten.BASE.VarParam.a     Mess.Messdaten.BASE.VarParam.b;
         Mess.Messdaten.O2.VarParam.a       Mess.Messdaten.O2.VarParam.b];

fprintf('Freie Parameter: %s\n', strjoin(namen(iFree), ', '));
fprintf('wsig = %s\n', mat2str(S.wsig));

%% 2. Auslegung des geplanten Experiments --------------------------------
VP.nx = 7;   VP.np = 10;
VP.p     = p;
VP.iFree = iFree;
VP.wsig  = S.wsig;              % identisch zum Fit

VP.T      = 30;                 % Horizont [h], wie die historischen Laeufe
VP.dt_off = 1.5;                % Abtastung offline [h]
VP.dt_on  = 0.1;                % Abtastung online (Base, DOT) [h]
VP.Vprobe = 0.03;               % Probenvolumen je Offline-Probe [L]

VP.DOTstern = max(TrainSet(2).O2.y);
VP.x0    = TrainSet(2).x0;
VP.x0(1) = 10;                  % V0 = 10 L  [ANNAHME -- Weight pruefen]
VP.x0(3) =  3 * VP.x0(1);       % cGlc0 = 3 g/L, knapp unter KS
VP.x0(5) = 0.5 * VP.x0(1);      % cPh0  = 0.5 g/L, damit Ph limitierend wird

% Feedkonzentrationen (Vorratsloesungen, nicht optimiert)
VP.cAm_in  = 30;   VP.cPh_in = 24;   VP.cGlc_in = 450;

% Zeitraster
VP.nu  = 10;                    % Abschnitte Glucose
VP.nuA = 3;                     % Abschnitte Ammonium
VP.nuP = 3;                     % Abschnitte Phosphat
VP.tu    = linspace(0, VP.T, VP.nu+1);
VP.t_off = (0:VP.dt_off:VP.T).';
VP.t_on  = (0:VP.dt_on :VP.T).';
VP.t_con = unique([(0:0.25:VP.T).'; VP.T]);   % Raster der Nebenbedingungen
VP.mapA  = ceil((1:VP.nu) * VP.nuA / VP.nu);  % grob -> fein
VP.mapP  = ceil((1:VP.nu) * VP.nuP / VP.nu);

VP.Probe.BatchAge = VP.t_off(2:end);
VP.Probe.Volumen  = VP.Vprobe * ones(numel(VP.t_off)-1, 1);

%% 3. Schranken -----------------------------------------------------------
% Die Feedmaxima sind an den Bedarf angepasst: bei ~200 g Biomasse braucht
% man ueber 30 h rund 0.6 L Ammonium- und 0.5 L Phosphatloesung.
VP.uGlc_max = 0.20;    % L/h
VP.uAm_max  = 0.05;    % L/h
VP.uPh_max  = 0.03;    % L/h

VP.Vmax     = 14.0;    % L  [ANNAHME: 15-L-Reaktor]
VP.cGlc_max = 8;  %35    % g/L, aus den historischen Laeufen
VP.DOTmin   = 20;      % %, keine O2-Limitierung (Modellgueltigkeit)
VP.cAm_min  = 0.10;    % g/L, deutlich ueber KAm
VP.cPh_min  = 0.02;    % g/L, immer noch 200x ueber KPh

%% 4. Ausgangslage --------------------------------------------------------
[relStd_alt, CV_alt, CVabs_alt, Fn_alt] = analyse_fim(FM_alt, p, iFree);
K_alt = crit_all(Fn_alt);
Corr_alt = corr_from_cov(CVabs_alt);

fprintf('\n--- VOR OVP (nur Altversuche) ---\n');
fprintf('%-9s %12s %10s\n','Param','Wert','rel. [%]');
for i = 1:nf
    fprintf('%-9s %12.4f %10.1f\n', namen{iFree(i)}, p(iFree(i)), 100*relStd_alt(i));
end
fprintf('A = %.4e | -logdet = %.3f | Kondition = %.3e\n', K_alt.A, K_alt.D, K_alt.CN);
print_worst_dir(Fn_alt, namen(iFree));

%% 5. Optimierung je Kriterium -------------------------------------------
z0 = [0.25*ones(1,VP.nu), ...
      (0.005/VP.uAm_max)*ones(1,VP.nuA), ...
      (0.002/VP.uPh_max)*ones(1,VP.nuP)];
nz = numel(z0);

opts = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 2000, ...
                    'FiniteDifferenceStepSize', 1e-3, ...
                    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8);

kritListe = {'A','D'};
Erg = struct('krit',{},'z',{},'relStd',{},'K',{});

for ic = 1:numel(kritListe)
    krit = kritListe{ic};
    fprintf('\n=== OVP mit %s-Kriterium ===\n', krit);
    tStart = tic;

    z_opt = fmincon(@(z) ovp_objective(z, VP, FM_alt, krit), z0, ...
                    [],[],[],[], zeros(1,nz), ones(1,nz), ...
                    @(z) ovp_constraints(z, VP), opts);

    [rS, ~, ~, Fn] = analyse_fim(FM_alt + fim_design(z_opt, VP), p, iFree);
    Erg(ic).krit   = krit;
    Erg(ic).z      = z_opt;
    Erg(ic).relStd = rS;
    Erg(ic).K      = crit_all(Fn);
    fprintf('Dauer %.1f min\n', toc(tStart)/60);
end

% Referenz: gleiches Experiment mit KONSTANTEN Feeds. Ohne diese Zeile
% misst der Vergleich nur den Wert eines zusaetzlichen Versuchs.
[relStd_ref, ~, ~, Fn_ref] = analyse_fim(FM_alt + fim_design(z0, VP), p, iFree);
K_ref = crit_all(Fn_ref);

%% 6. Ergebnistabellen ----------------------------------------------------
fprintf('\n=== KRITERIEN (Zeile = optimiertes Kriterium) ===\n');
fprintf('%-12s %14s %14s %14s\n','Entwurf','A','-logdet','Kondition');
fprintf('%-12s %14.4e %14.3f %14.3e\n','nur alt',  K_alt.A, K_alt.D, K_alt.CN);
fprintf('%-12s %14.4e %14.3f %14.3e\n','konstant', K_ref.A, K_ref.D, K_ref.CN);
for ic = 1:numel(Erg)
    fprintf('%-12s %14.4e %14.3f %14.3e\n', Erg(ic).krit, ...
            Erg(ic).K.A, Erg(ic).K.D, Erg(ic).K.CN);
end

fprintf('\n=== RELATIVE STANDARDABWEICHUNGEN [%%] ===\n');
fprintf('%-12s', 'Entwurf'); fprintf('%10s', namen{iFree}); fprintf('\n');
fprintf('%-12s', 'nur alt');  fprintf('%10.1f', 100*relStd_alt); fprintf('\n');
fprintf('%-12s', 'konstant'); fprintf('%10.1f', 100*relStd_ref); fprintf('\n');
for ic = 1:numel(Erg)
    fprintf('%-12s', Erg(ic).krit);
    fprintf('%10.1f', 100*Erg(ic).relStd); fprintf('\n');
end

%% 7. Bester Entwurf ------------------------------------------------------
kritWahl = 'A'; %'D';
icBest   = find(strcmp({Erg.krit}, kritWahl), 1);
z_best   = Erg(icBest).z;
[FM_neu, X, t_sim] = fim_design(z_best, VP);
[~, ~, CVabs_neu, Fn_neu] = analyse_fim(FM_alt + FM_neu, p, iFree);
Corr_neu = corr_from_cov(CVabs_neu);

fprintf('\n--- Schlechteste Richtung nach OVP (%s) ---\n', kritWahl);
print_worst_dir(Fn_neu, namen(iFree));

[uG, uA, uP] = unpack_u(z_best, VP);
fprintf('\n--- Optimale Feeds (%s) ---\n', kritWahl);
fprintf('u_Glc [L/h]: %s\n', mat2str(uG, 3));
fprintf('u_Am  [L/h]: %s\n', mat2str(VP.uAm_max*z_best(VP.nu+1:VP.nu+VP.nuA), 3));
fprintf('u_Ph  [L/h]: %s\n', mat2str(VP.uPh_max*z_best(VP.nu+VP.nuA+1:end), 3));
fprintf('An der Schranke: %d von %d Variablen\n', ...
        nnz(z_best > 0.999 | z_best < 0.001), nz);

save(fullfile(DATEN,'FIM','OVP_Modell3.mat'), ...
     'Erg','K_alt','K_ref','relStd_alt','relStd_ref','VP','iFree','z_best');

%% 8. Abbildungen ---------------------------------------------------------
FS = 16;  FSA = 12;
V = X(:,1);

figure('Name','OVP Modell 3 -- optimale Feeds','Position',[100 40 950 800]);

% Alle drei Feeds in einem Plot. Am/Ph liegen eine Groessenordnung unter
% Glucose, deshalb auf der rechten Achse.
subplot(4,1,1);
yyaxis left
stairs(VP.tu, [uG uG(end)], 'LineWidth', 2);
ylabel('u_{Glc} [L/h]','FontSize',FS);
yyaxis right
stairs(VP.tu, [uA uA(end)], 'LineWidth', 2); hold on;
stairs(VP.tu, [uP uP(end)], '-.', 'LineWidth', 2);
ylabel('u_{Am}, u_{Ph} [L/h]','FontSize',FS);
legend('u_{Glc}','u_{Am}','u_{Ph}','Location','best');
set(gca,'FontSize',FSA); grid on;
title(sprintf('Optimale Feeds (%s-Kriterium)', kritWahl), 'FontSize', FS);

subplot(4,1,2);
plot(t_sim, X(:,3)./V, 'LineWidth', 2);
ylabel('c_{Glc} [g/L]','FontSize',FS); set(gca,'FontSize',FSA); grid on;

% Am/Ph liegen an ihren unteren Schranken -- gegen cX unsichtbar,
% deshalb eigene Achse rechts.
subplot(4,1,3);
yyaxis left
plot(t_sim, X(:,2)./V, 'LineWidth', 2);
ylabel('c_X [g/L]','FontSize',FS);
yyaxis right
plot(t_sim, X(:,4)./V, 'LineWidth', 1.5); hold on;
plot(t_sim, X(:,5)./V, '-.', 'LineWidth', 1.5);
yline(VP.cAm_min, ':k'); yline(VP.cPh_min, ':k');
ylabel('c_{Am}, c_{Ph} [g/L]','FontSize',FS);
ylim([0, 1.3*max([VP.cAm_min, VP.cPh_min, ...
                  max(X(:,4)./V), max(X(:,5)./V)])]);
legend('c_X','c_{Am}','c_{Ph}','Location','best');
set(gca,'FontSize',FSA); grid on;

subplot(4,1,4);
yyaxis left;  plot(t_sim, V, 'LineWidth', 2);      ylabel('V [L]','FontSize',FS);
yyaxis right; plot(t_sim, X(:,7), 'LineWidth', 2); ylabel('DOT [%]','FontSize',FS);
xlabel('Zeit [h]','FontSize',FS); set(gca,'FontSize',FSA); grid on;

build_heatmap(Corr_alt, namen(iFree), 'Korrelation vor OVP');
build_heatmap(Corr_neu, namen(iFree), sprintf('Korrelation nach OVP (%s)', kritWahl));

figure('Name','Feedprofile A vs D','Position',[100 100 900 600]);
for ic = 1:numel(Erg)
    [gA, aA, pA] = unpack_u(Erg(ic).z, VP);
    subplot(3,1,1); stairs(VP.tu, [gA gA(end)], 'LineWidth', 2); hold on;
    subplot(3,1,2); stairs(VP.tu, [aA aA(end)], 'LineWidth', 2); hold on;
    subplot(3,1,3); stairs(VP.tu, [pA pA(end)], 'LineWidth', 2); hold on;
end
subplot(3,1,1); ylabel('u_{Glc} [L/h]','FontSize',FS);
legend({Erg.krit},'Location','best'); grid on; set(gca,'FontSize',FSA);
subplot(3,1,2); ylabel('u_{Am} [L/h]','FontSize',FS); grid on; set(gca,'FontSize',FSA);
subplot(3,1,3); ylabel('u_{Ph} [L/h]','FontSize',FS); xlabel('Zeit [h]','FontSize',FS);
grid on; set(gca,'FontSize',FSA);

%% 9. Alles als SVG speichern --------------------------------------------
bildordner = fullfile(projectRoot,'Bilder','OVP');
figs = findobj('Type','figure');
[~, ord] = sort([figs.Number]);
for k = ord(:).'
    save_fig(figs(k), figs(k).Name, bildordner);
end


%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function J = ovp_objective(z, VP, FM_alt, krit)
    try
        FM_neu = fim_design(z, VP);
    catch
        J = 1e12; return;
    end
    iF = VP.iFree;  pf = VP.p(iF);
    Fn = diag(pf) * (FM_alt(iF,iF) + FM_neu(iF,iF)) * diag(pf);
    if any(~isfinite(Fn(:))) || rcond(Fn) < 1e-14, J = 1e12; return; end
    K = crit_all(Fn);
    switch krit
        case 'A', J = K.A;
        case 'D', J = K.D;
        otherwise, error('Unbekanntes Kriterium %s', krit);
    end
    if ~isfinite(J), J = 1e12; end
end


function K = crit_all(Fn)
% A = trace(inv Fn), D = -logdet(Fn). logdet ueber Cholesky, nicht ueber
% det: die Eigenwerte spannen mehrere Groessenordnungen.
    ev = max(sort(eig((Fn+Fn.')/2),'descend'), eps);
    [R, flag] = chol((Fn+Fn.')/2);
    if flag == 0, K.D = -2*sum(log(diag(R))); else, K.D = -sum(log(ev)); end
    K.A  = sum(1./ev);
    K.CN = ev(1)/ev(end);
end


function [c, ceq] = ovp_constraints(z, VP)
% Zustandsbeschraenkungen als VEKTOR ueber der Zeit. Ein max ueber die
% Trajektorie waere nicht differenzierbar -- SQP wuerde die Restriktion
% erst sehen, wenn sie schon verletzt ist.
    ceq = [];
    try
        X = cached_states(z, VP);
    catch
        c = 1e3 * ones(5*numel(VP.t_con), 1); return;
    end
    V = X(:,1);
    c = [  V         - VP.Vmax;
           X(:,3)./V - VP.cGlc_max;
           VP.DOTmin - X(:,7);
           VP.cAm_min - X(:,4)./V;
           VP.cPh_min - X(:,5)./V ];
end


function X = cached_states(z, VP)
% fmincon ruft Ziel- und Nebenbedingungsfunktion am selben Punkt auf.
    persistent zLast Xlast
    if ~isempty(zLast) && isequal(zLast, z), X = Xlast; return; end
    X = sim_core(VP.t_con, VP.x0(:), build_u(z,VP), VP.p, VP.DOTstern, ...
                 VP.Probe, @Modell3_woEtOH_10p, VP.nx, 0);
    zLast = z;  Xlast = X;
end


function [FM, X, t_sim] = fim_design(z, VP)
% FIM des geplanten Experiments. Gewichtung wie im Guetefunktional:
% wsig(i) je Kanal und Mittelung ueber die Punktzahl (Faktor 1/n).
    nx = VP.nx;  np = VP.np;
    t_sim = unique([VP.t_off; VP.t_on; VP.tu(:)]);
    X_ext = sim_core(t_sim, [VP.x0(:); zeros(nx*np,1)], build_u(z,VP), VP.p, ...
                     VP.DOTstern, VP.Probe, @Modell3_woEtOH_XP, nx, np);
    X = X_ext(:, 1:nx);

    % {Index in der Messgleichung, Zeitraster, Zeile in VP.ab}
    kan = { 1, VP.t_off, 1;  2, VP.t_off, 2;  3, VP.t_off, 3; ...
            4, VP.t_off, 4;  5, VP.t_on,  5;  6, VP.t_on,  6 };

    FM = zeros(np, np);
    for i = 1:6
        if VP.wsig(i) == 0, continue; end
        iOut = kan{i,1};
        [tf, idx] = ismember(kan{i,2}, t_sim);
        if any(~tf), error('Zeitpunkt nicht im Simulationsraster.'); end

        n  = numel(idx);
        Sy = zeros(n, np);   yv = zeros(n, 1);
        for k = 1:n
            xk   = X_ext(idx(k), 1:nx).';
            XP_k = reshape(X_ext(idx(k), nx+1:end), nx, np);
            dhdx    = Modell3_dmgldx(xk);
            Sy(k,:) = dhdx(iOut,:) * XP_k;
            yk = Modell3_mgl(xk);   yv(k) = yk(iOut);
        end

        a = VP.ab(kan{i,3},1);   b = VP.ab(kan{i,3},2);
        if iOut == 5
            v = a.*yv/1000 + b/1e6;      % Base: Doku mL -> Modell L
        else
            v = a.*yv + b;               % O2 braucht keine Umrechnung
        end
        v = max(v, eps);

        for k = 1:n
            FM = FM + (VP.wsig(i)/n) * (Sy(k,:).' * Sy(k,:)) / v(k);
        end
    end
end


function [uG, uA, uP] = unpack_u(z, VP)
    nG = VP.nu;  nA = VP.nuA;
    uG = VP.uGlc_max * reshape(z(1:nG), 1, []);
    uA = VP.uAm_max  * reshape(z(nG + VP.mapA), 1, []);
    uP = VP.uPh_max  * reshape(z(nG + nA + VP.mapP), 1, []);
end


function u = build_u(z, VP)
    [uG, uA, uP] = unpack_u(z, VP);
    u = zeros(10, numel(VP.tu));
    u(1,:) = VP.tu;
    u(2,:) = [uA uA(end)];   u(3,:) = VP.cAm_in;
    u(4,:) = [uP uP(end)];   u(5,:) = VP.cPh_in;
    u(6,:) = [uG uG(end)];   u(7,:) = VP.cGlc_in;
end


function X = sim_core(t, x0, u, p, DOTstern, Probe, rhs, nx, np)
% Segmentweise Integration mit Probenahme-Spruengen. NonNegative nur im
% reinen Zustandslauf -- im Sensitivitaetslauf waere die Klammerung nicht
% differenzierbar, dort halten die Nebenbedingungen die Zustaende positiv.
    t  = t(:);
    tp = Probe.BatchAge(:);  vp = Probe.Volumen(:);
    tf = u(1,:).';           t0 = u(1,1);
    I  = tp >= t0 & tp <= max(t);  tp = tp(I);  vp = vp(I);
    tf = tf(tf >= t0 & tf <= max(t));
    tb = unique([t0; tf; tp; max(t)]);

    x = x0(:);
    X = zeros(numel(t), numel(x0));
    if np == 0
        opt = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5,'NonNegative',1:nx);
    else
        opt = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5);
    end

    X(t == t0,:) = repmat(x.', nnz(t == t0), 1);
    for j = find(tp == t0).', x = probe_any(x, vp(j), nx, np); end

    for k = 2:numel(tb)
        ta = tb(k-1);  te = tb(k);
        sel  = t > ta & t <= te;
        tout = unique([ta; t(sel); te]);
        [~, Xs] = ode15s(@(tt,xx) rhs(tt,xx,u,p,DOTstern), tout, x, opt);
        if numel(tout) == 2, Xs = Xs([1 end],:); end
        if size(Xs,1) ~= numel(tout), error('Integration unvollstaendig.'); end
        idxT = find(sel);
        if ~isempty(idxT)
            [~, loc] = ismember(t(idxT), tout);
            X(idxT,:) = Xs(loc,:);
        end
        x = Xs(end,:).';
        for j = find(tp == te).', x = probe_any(x, vp(j), nx, np); end
    end
end


function xe = probe_any(xe, Vp, nx, np)
% Probenahme auf Zustand und Sensitivitaeten: x+ = g(x), also XP+ = dg/dx*XP.
    x = xe(1:nx);
    if np == 0, xe = probe_m3(x, Vp); return; end
    XP = reshape(xe(nx+1:end), nx, np);
    Jg = zeros(nx);
    for j = 1:nx
        h  = 1e-6 * max(abs(x(j)), 1);
        xp = x; xp(j) = xp(j) + h;
        xm = x; xm(j) = xm(j) - h;
        Jg(:,j) = (probe_m3(xp,Vp) - probe_m3(xm,Vp)) / (2*h);
    end
    xe = [reshape(probe_m3(x,Vp), [], 1); reshape(Jg*XP, [], 1)];
end


function [relStd, CV, CVabs, Fn] = analyse_fim(FM, p, iFree)
% Durch die Normierung sind die Diagonalelemente von CV bereits relative
% Varianzen; CVabs ist die Kovarianz in physikalischen Einheiten.
    pf = p(iFree);
    Fn = diag(pf) * FM(iFree,iFree) * diag(pf);
    if rcond(Fn) < 1e-14, CV = pinv(Fn); else, CV = Fn \ eye(numel(iFree)); end
    relStd = sqrt(abs(diag(CV))).';
    CVabs  = diag(pf) * CV * diag(pf);
end


function C = corr_from_cov(CV)
    sd = sqrt(abs(diag(CV)));
    C  = diag(1./max(sd,eps)) * CV * diag(1./max(sd,eps));
end


function print_worst_dir(Fn, nam)
% Eigenvektor zum kleinsten Eigenwert: welche Parameterkombination ist am
% schlechtesten bestimmt?
    [EV, EW] = eig((Fn+Fn.')/2);
    [~, imin] = min(diag(EW));
    v = EV(:,imin);
    [~, ord] = sort(abs(v), 'descend');
    fprintf('Schlechteste Richtung: ');
    for k = 1:min(3, numel(v)), fprintf('%+.2f*%s ', v(ord(k)), nam{ord(k)}); end
    fprintf('\n');
end


function build_heatmap(corr, params, titel)
n = 256;
blue = [0.0000 0.4470 0.7410];  white = [1 1 1];  red = [0.8500 0.3250 0.0980];
cmap = [ [linspace(blue(1),white(1),n/2)' linspace(blue(2),white(2),n/2)' ...
          linspace(blue(3),white(3),n/2)'];
         [linspace(white(1),red(1),n/2)'  linspace(white(2),red(2),n/2)'  ...
          linspace(white(3),red(3),n/2)'] ];

figure('Color','w','Name',titel);
h = heatmap(params, params, corr);
h.Title = titel;
h.XLabel = 'Parameter';   h.YLabel = 'Parameter';
h.ColorLimits = [-1 1];   h.Colormap = cmap;
h.CellLabelFormat = '%.2f';   h.FontSize = 11;
end


function save_fig(fig, name, ordner)
% Speichert eine Figure als SVG.
    if ~exist(ordner,'dir'), mkdir(ordner); end
    name = regexprep(strtrim(name), '[^\w\-]', '_');
    if isempty(name), name = sprintf('Figure_%d', fig.Number); end
    ziel = fullfile(ordner, [name '.svg']);

    set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
    drawnow;
    try
        exportgraphics(fig, ziel, 'ContentType', 'vector');
    catch
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, ziel, '-dsvg');
    end
    fprintf('  gespeichert: %s\n', ziel);
end