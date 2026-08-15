% main_06_ovp_multikriterium.m
%
% Optimale Versuchsplanung (OVP) fuer Modell 3 -- erweiterte Fassung.
%
% Aenderungen gegenueber main_06_ovp_neu.m:
%   1. Optimiert werden DREI Feeds (Glucose, Ammonium, Phosphat), nicht nur
%      u_Glc. Am/Ph laufen auf einem groeberen Zeitraster als Glc, weil ihr
%      informativer Gehalt ein Umschalten ist und keine feine Form.
%   2. Entscheidungsvariablen sind normiert auf [0,1]. Sonst ist ein
%      fmincon-Differenzenschritt fuer u_Ph relativ 100x groesser als fuer
%      u_Glc und die Gradienten sind nicht vergleichbar.
%   3. Nebenbedingungen als VEKTOR ueber der Zeit statt max/min. Ein max
%      ueber eine Trajektorie ist nicht differenzierbar -- der FD-Gradient
%      ist fast ueberall null und SQP "sieht" die Restriktion nicht.
%   4. Untere Schranken fuer cAm und cPh. Das Modell klammert nicht, ohne
%      diese Schranken kann der Optimierer die Kultur aushungern und die
%      Monod-Terme kippen das Vorzeichen.
%   5. Vier Kriterien (A, D, E, mod-E) statt nur A, jeweils optimiert und
%      anschliessend alle vier auf allen Loesungen ausgewertet (Kreuztabelle).
%   6. Optional sequentielle Planung mehrerer Experimente (Cramer-Rao ueber
%      die Summe der FIM, vgl. Herold et al. 2017, Gl. 17).
%
% Hinweis zur Normierung: A, E und mod-E brauchen die diag(p)-Skalierung,
% sonst dominiert YXO (~246) gegenueber YB_Am (~0.02). D ist invariant
% gegenueber linearer Reparametrisierung -- die Skalierung verschiebt das
% D-Kriterium nur um eine Konstante und aendert das Optimum nicht.

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fullfile(projectRoot,'..', '..');
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

FM_alt = F0.FM;                       % 10x10 aus den vier Trainingslaeufen

% FM_alt und p muessen aus DEMSELBEN Fit stammen
if max(abs(F0.p_opt(:) - p)) > 1e-9
    error(['FIM_Modell3_multi.mat wurde mit anderen Parametern gerechnet. ', ...
           'Erst main_05_uncertainty_FIM_m3_multi.m neu laufen lassen.']);
end

% Varianzparameter a,b (gleiche Quelle und Umrechnung wie im Preprocessing)
Mess = load(fullfile(DATEN,'MessDaten_SPI1_Projekt','Mess_RamScDef10.mat')).Mess;
VP.ab = [Mess.Messdaten.Biomasse.VarParam.a Mess.Messdaten.Biomasse.VarParam.b;
         Mess.Messdaten.Glucose.VarParam.a  Mess.Messdaten.Glucose.VarParam.b;
         Mess.Messdaten.Ammonium.VarParam.a Mess.Messdaten.Ammonium.VarParam.b;
         Mess.Messdaten.Phosphat.VarParam.a Mess.Messdaten.Phosphat.VarParam.b;
         Mess.Messdaten.BASE.VarParam.a     Mess.Messdaten.BASE.VarParam.b;
         Mess.Messdaten.O2.VarParam.a       Mess.Messdaten.O2.VarParam.b];

fprintf('Freie Parameter: %s\n', strjoin(namen(iFree), ', '));

%% 2. Auslegung des geplanten Experiments --------------------------------
VP.nx = 7;   VP.np = 10;
VP.p        = p;
VP.iFree    = iFree;
VP.T        = 15;          % Horizont [h]
VP.dt_off   = 1.5;         % Abtastung offline (Bio, Glc, Am, Ph) [h]
VP.dt_on    = 0.1;         % Abtastung online (Base, DOT) [h]
VP.dt_base  = 2.0;         % Ausduennung der Base-Inkremente (wie im Fit)
VP.Vprobe   = 0.03;        % Probenvolumen je Offline-Probe [L]
VP.DOTstern = max(TrainSet(2).O2.y);
VP.x0       = TrainSet(2).x0;
VP.x0(1)    = 10;                    % V0 = 10 L
VP.x0(3)    = 3 * VP.x0(1);          % cGlc0 ~ 3 g/L, sonst laeuft der halbe
                                     % Horizont fuer den Abbau der Startcharge

% Feedkonzentrationen (Vorratsloesungen, nicht optimiert)
VP.cAm_in  = 30;
VP.cPh_in  = 24;
VP.cGlc_in = 450;

% Zeitraster
VP.nu  = 10;               % Abschnitte Glucose-Feed
VP.nuA = 3;                % Abschnitte Ammonium-Feed (grober)
VP.nuP = 3;                % Abschnitte Phosphat-Feed (grober)
VP.tu    = linspace(0, VP.T, VP.nu+1);
VP.t_off = (0:VP.dt_off:VP.T).';
VP.t_on  = (0:VP.dt_on :VP.T).';
VP.t_con = unique([(0:0.25:VP.T).'; VP.T]);   % Raster der Nebenbedingungen

% Abbildung grobes -> feines Feedraster (Nullter-Ordnung-Halt)
VP.mapA = ceil((1:VP.nu) * VP.nuA / VP.nu);
VP.mapP = ceil((1:VP.nu) * VP.nuP / VP.nu);

% Probenahmen an den Offline-Zeitpunkten (ausser t = 0)
VP.Probe.BatchAge = VP.t_off(2:end);
VP.Probe.Volumen  = VP.Vprobe * ones(numel(VP.t_off)-1, 1);

%% 3. Schranken und Nebenbedingungen -------------------------------------
VP.uGlc_max = 0.20;    % L/h
VP.uAm_max  = 0.02;    % L/h
VP.uPh_max  = 0.01;    % L/h

VP.Vmax     = 14.0;    % L, Reaktorgrenze
VP.cGlc_max = 35;      % g/L, aus vergangenen Versuchen
VP.DOTmin   = 20;      % %, keine O2-Limitierung erzwingen

% Untere Schranken fuer Am/Ph. Bewusst deutlich ueber KAm/KPh: sonst
% entkoppelt die Planung YXS/YAmX zwar besser, haengt dafuer aber an einem
% fixierten, nicht identifizierbaren Halbsaettigungswert. Wer die
% aggressive Variante will, setzt hier ~2*KAm ein.
VP.cAm_min = 0.10;     % g/L
VP.cPh_min = 0.05;     % g/L

%% 4. Ausgangslage: was wissen wir ohne neues Experiment? ----------------
[relStd_alt, CN_alt, CV_alt, CVabs_alt, Fn_alt] = analyse_fim(FM_alt, p, iFree);
K_alt = crit_all(Fn_alt);

fprintf('\n--- VOR OVP (nur Altversuche) ---\n');
fprintf('%-9s %12s %10s\n','Param','Wert','rel. [%]');
for i = 1:nf
    fprintf('%-9s %12.4f %10.1f\n', namen{iFree(i)}, p(iFree(i)), 100*relStd_alt(i));
end
fprintf('A = %.4e | -logdet = %.3f | 1/lmin = %.4e | Kondition = %.3e\n', ...
        K_alt.A, K_alt.D, K_alt.E, K_alt.ME);
print_worst_dir(Fn_alt, namen(iFree));

%% 5. Startpunkt und Optimiereroptionen ----------------------------------
% z ist normiert: z = [zGlc(1..nu) zAm(1..nuA) zPh(1..nuP)] in [0,1]
z0 = [0.25*ones(1,VP.nu), ...
      (0.005/VP.uAm_max)*ones(1,VP.nuA), ...
      (0.002/VP.uPh_max)*ones(1,VP.nuP)];
nz  = numel(z0);
zLB = zeros(1,nz);
zUB = ones(1,nz);

opts = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                    'MaxFunctionEvaluations', 2000, ...
                    'FiniteDifferenceStepSize', 1e-3, ...   % z ist O(1)
                    'OptimalityTolerance', 1e-6, ...
                    'StepTolerance', 1e-8, ...
                    'UseParallel', false);   % ggf. auf true setzen

%% 6. Je Kriterium eine Optimierung --------------------------------------
kritListe = {'A','D','E','ME'};
kritText  = {'A (Summe der Varianzen)', 'D (Volumen des Ellipsoids)', ...
             'E (schlechteste Richtung)', 'mod-E (Kondition)'};

Erg = struct('krit',{},'z',{},'FM_neu',{},'relStd',{},'K',{},'t',{});

for ic = 1:numel(kritListe)
    krit = kritListe{ic};
    fprintf('\n=== OVP mit %s-Kriterium: %s ===\n', krit, kritText{ic});
    tStart = tic;

    obj = @(z) ovp_objective(z, VP, FM_alt, krit);
    con = @(z) ovp_constraints(z, VP);
    z_opt = fmincon(obj, z0, [], [], [], [], zLB, zUB, con, opts);

    FM_neu = fim_design(z_opt, VP);
    [rS, ~, ~, ~, Fn] = analyse_fim(FM_alt + FM_neu, p, iFree);

    Erg(ic).krit   = krit;
    Erg(ic).z      = z_opt;
    Erg(ic).FM_neu = FM_neu;
    Erg(ic).relStd = rS;
    Erg(ic).K      = crit_all(Fn);
    Erg(ic).t      = toc(tStart);
    fprintf('Dauer %.1f min\n', Erg(ic).t/60);
end

%% 7. Kreuztabelle: jedes Kriterium auf jeder Loesung --------------------
% Der interessante Befund ist die Tabelle, nicht eine einzelne Zahl. Bei
% schlecht konditionierter FIM gilt A ~ 1/lmin, dann fallen A- und
% E-optimaler Entwurf praktisch zusammen.
fprintf('\n=== KREUZTABELLE (Zeile = optimiertes Kriterium) ===\n');
fprintf('%-12s %12s %12s %12s %12s\n','Entwurf','A','-logdet','1/lmin','Kondition');
fprintf('%-12s %12.4e %12.3f %12.4e %12.3e\n','nur alt', ...
        K_alt.A, K_alt.D, K_alt.E, K_alt.ME);
for ic = 1:numel(Erg)
    fprintf('%-12s %12.4e %12.3f %12.4e %12.3e\n', Erg(ic).krit, ...
            Erg(ic).K.A, Erg(ic).K.D, Erg(ic).K.E, Erg(ic).K.ME);
end

fprintf('\n=== RELATIVE STANDARDABWEICHUNGEN [%%] ===\n');
fprintf('%-12s', 'Entwurf'); fprintf('%10s', namen{iFree}); fprintf('\n');
fprintf('%-12s', 'nur alt');  fprintf('%10.1f', 100*relStd_alt); fprintf('\n');
for ic = 1:numel(Erg)
    fprintf('%-12s', Erg(ic).krit);
    fprintf('%10.1f', 100*Erg(ic).relStd); fprintf('\n');
end

% Referenzentwurf: gleiches Experiment, aber NICHT optimierte Feeds.
% Ohne diese Zeile misst der Vergleich nur den Wert eines ZUSAETZLICHEN
% Versuchs, nicht den Wert der Optimierung.
FM_ref = fim_design(z0, VP);
[relStd_ref, ~, ~, ~, Fn_ref] = analyse_fim(FM_alt + FM_ref, p, iFree);
K_ref = crit_all(Fn_ref);
fprintf('%-12s', 'konstant');  fprintf('%10.1f', 100*relStd_ref); fprintf('\n');
fprintf('\nReferenz (konstante Feeds): A = %.4e | -logdet = %.3f\n', K_ref.A, K_ref.D);

%% 8. Bester Entwurf: Kennwerte und schlechteste Richtung ----------------
kritWahl = 'D';                       % Entwurf fuer Plots und Sequenz
icBest   = find(strcmp({Erg.krit}, kritWahl), 1);
z_best   = Erg(icBest).z;
[FM_best, X, t_sim] = fim_design(z_best, VP);
[~, ~, CV_neu, CVabs_neu, Fn_neu] = analyse_fim(FM_alt + FM_best, p, iFree);

fprintf('\n--- Schlechteste Richtung nach OVP (%s) ---\n', kritWahl);
print_worst_dir(Fn_neu, namen(iFree));

pf = p(iFree);
[Corr_neu, ~, ~, ~, ~, ~] = Parameteranalyse(CVabs_neu, pf);
fprintf('\n--- Korrelationsmatrix nach OVP ---\n');
fprintf('%9s',''); fprintf('%9s', namen{iFree}); fprintf('\n');
for i = 1:nf
    fprintf('%9s', namen{iFree(i)});
    fprintf('%9.2f', Corr_neu(i,:));  fprintf('\n');
end

%% 9. Sequentielle Planung mehrerer Experimente --------------------------
% Wie in Herold et al. 2017, Tab. 2: jedes weitere optimal geplante
% Experiment senkt die Unsicherheit. Die FIM addieren sich (Gl. 17), also
% wird jeder Folgeentwurf gegen die bereits akkumulierte Information
% optimiert. Beantwortet die Frage "wie viele Versuche braeuchte man?"
doSequenz = false;                    % teuer: nPlan x eine volle Optimierung
nPlan     = 3;

if doSequenz
    FM_acc  = FM_alt;
    seqStd  = zeros(nPlan, nf);
    seqZ    = zeros(nPlan, nz);
    for k = 1:nPlan
        fprintf('\n=== Sequentielle Planung, Experiment %d ===\n', k);
        obj = @(z) ovp_objective(z, VP, FM_acc, kritWahl);
        con = @(z) ovp_constraints(z, VP);
        zk  = fmincon(obj, z0, [], [], [], [], zLB, zUB, con, opts);
        FM_acc = FM_acc + fim_design(zk, VP);
        [seqStd(k,:), ~, ~, ~, ~] = analyse_fim(FM_acc, p, iFree);
        seqZ(k,:) = zk;
    end
    fprintf('\n=== SEQUENZ: relative Standardabweichung [%%] ===\n');
    fprintf('%-12s', 'Stand'); fprintf('%10s', namen{iFree}); fprintf('\n');
    fprintf('%-12s', 'nur alt'); fprintf('%10.1f', 100*relStd_alt); fprintf('\n');
    for k = 1:nPlan
        fprintf('%-12s', sprintf('+ Exp %d', k));
        fprintf('%10.1f', 100*seqStd(k,:)); fprintf('\n');
    end
else
    seqStd = []; seqZ = [];
end

%% 10. Speichern ----------------------------------------------------------
if ~exist(fullfile(DATEN,'FIM'),'dir'), mkdir(fullfile(DATEN,'FIM')); end
save(fullfile(DATEN,'FIM','OVP_Modell3_multikriterium.mat'), ...
     'Erg','K_alt','K_ref','relStd_alt','relStd_ref','VP','iFree', ...
     'z_best','kritWahl','seqStd','seqZ');

%% 11. Abbildungen --------------------------------------------------------
[uG, uA, uP] = unpack_u(z_best, VP);
V  = X(:,1);
FS = 16; FSA = 12;

figure('Name','OVP Modell 3 -- optimale Feeds','Position',[100 40 950 950]);
subplot(5,1,1);
stairs(VP.tu, [uG uG(end)], 'LineWidth', 2);
ylabel('u_{Glc} [L/h]','FontSize',FS); set(gca,'FontSize',FSA); grid on;
title(sprintf('Optimale Feeds (%s-Kriterium)', kritWahl), 'FontSize', FS);

subplot(5,1,2);
stairs(VP.tu, [uA uA(end)], 'LineWidth', 2); hold on;
stairs(VP.tu, [uP uP(end)], 'LineWidth', 2);
ylabel('u [L/h]','FontSize',FS); legend('u_{Am}','u_{Ph}','Location','best');
set(gca,'FontSize',FSA); grid on;

subplot(5,1,3);
plot(t_sim, X(:,3)./V, 'LineWidth', 2); hold on;
yline(p(2),'--r','K_S');
ylabel('c_{Glc} [g/L]','FontSize',FS); set(gca,'FontSize',FSA); grid on;

subplot(5,1,4);
plot(t_sim, X(:,2)./V, 'LineWidth', 2); hold on;
plot(t_sim, X(:,4)./V, 'LineWidth', 1.5);
plot(t_sim, X(:,5)./V, 'LineWidth', 1.5);
yline(VP.cAm_min,'--k');
ylabel('c [g/L]','FontSize',FS); legend('c_X','c_{Am}','c_{Ph}','Location','best');
set(gca,'FontSize',FSA); grid on;

subplot(5,1,5);
yyaxis left;  plot(t_sim, V, 'LineWidth', 2);      ylabel('V [L]','FontSize',FS);
yyaxis right; plot(t_sim, X(:,7), 'LineWidth', 2); ylabel('DOT [%]','FontSize',FS);
xlabel('Zeit [h]','FontSize',FS); set(gca,'FontSize',FSA); grid on;

figure('Name','Parameterunsicherheit je Kriterium');
B = 100*[relStd_alt(:), relStd_ref(:), cell2mat(cellfun(@(r) r(:), {Erg.relStd}, ...
         'UniformOutput', false))];
bar(B);
set(gca,'XTickLabel', namen(iFree), 'FontSize', FSA);
ylabel('relative Standardabweichung [%]','FontSize',FS);
legend([{'vor OVP','konstant'}, {Erg.krit}], 'Location','best');
title('Wirkung der Kriterien','FontSize',FS); grid on;

% Unsicherheitsellipse fuer das staerkst korrelierte Paar YXS / YAmX
idx2 = [3 4];
figure('Name','Unsicherheitsellipse vor/nach OVP');
h1 = plot_gaussian_ellipsoid(pf(idx2), CVabs_alt(idx2,idx2), 1); hold on;
h2 = plot_gaussian_ellipsoid(pf(idx2), CVabs_neu(idx2,idx2), 1);
set(h1,'Color',[0.30 0.30 0.30],'LineWidth',1.5);
set(h2,'Color',[0.85 0.33 0.10],'LineWidth',2.0);
xlabel('Y_{XS}','FontSize',FS); ylabel('Y_{AmX}','FontSize',FS);
legend([h1 h2],'vor OVP',['nach OVP (' kritWahl ')'],'Location','best');
title('1\sigma-Unsicherheitsellipse','FontSize',FS);
set(gca,'FontSize',FSA); grid on;


%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function J = ovp_objective(z, VP, FM_alt, krit)
% Zielfunktion: gewaehltes Kriterium auf der normierten Gesamt-FIM.
    try
        FM_neu = fim_design(z, VP);
    catch
        J = 1e12; return;                 % Integration gescheitert
    end
    iF = VP.iFree;
    pf = VP.p(iF);
    Fn = diag(pf) * (FM_alt(iF,iF) + FM_neu(iF,iF)) * diag(pf);
    if any(~isfinite(Fn(:))) || rcond(Fn) < 1e-14
        J = 1e12; return;
    end
    K = crit_all(Fn);
    switch krit
        case 'A',  J = K.A;
        case 'D',  J = K.D;               % = -logdet(Fn)
        case 'E',  J = K.E;               % = 1/lambda_min
        case 'ME', J = K.ME;              % = lambda_max/lambda_min
        otherwise, error('Unbekanntes Kriterium %s', krit);
    end
    if ~isfinite(J), J = 1e12; end
end


function K = crit_all(Fn)
% Alle vier Kriterien auf einer normierten FIM. Alle als MINIMIERUNGSziel.
% logdet ueber Cholesky, nicht ueber det: die Eigenwerte spannen mehrere
% Groessenordnungen, det laeuft ueber bzw. unter.
    ev = sort(eig((Fn+Fn.')/2), 'descend');
    ev = max(ev, eps);
    [R, flag] = chol((Fn+Fn.')/2);
    if flag == 0
        logdetF = 2*sum(log(diag(R)));
    else
        logdetF = sum(log(ev));           % Notfall bei numerisch indefinit
    end
    K.A  = sum(1./ev);                    % = trace(inv(Fn))
    K.D  = -logdetF;                      % D: Volumen des Ellipsoids
    K.E  = 1/ev(end);                     % E: schlechteste Richtung
    K.ME = ev(1)/ev(end);                 % mod-E: Kondition, Form statt Groesse
    K.ev = ev;
end


function [c, ceq] = ovp_constraints(z, VP)
% Zustandsbeschraenkungen als VEKTOR ueber der Zeit (nicht als max/min).
% Nur Zustaende, ohne Sensitivitaeten -- deutlich billiger als die FIM.
    ceq = [];
    try
        X = cached_states(z, VP);
    catch
        c = 1e3 * ones(5*numel(VP.t_con), 1); return;
    end
    V = X(:,1);
    c = [  V            - VP.Vmax;          % V    <= Vmax
           X(:,3)./V    - VP.cGlc_max;      % cGlc <= cGlc_max
           VP.DOTmin    - X(:,7);           % DOT  >= DOTmin
           VP.cAm_min   - X(:,4)./V;        % cAm  >= cAm_min
           VP.cPh_min   - X(:,5)./V ];      % cPh  >= cPh_min
end


function X = cached_states(z, VP)
% fmincon ruft Ziel- und Nebenbedingungsfunktion am selben Punkt auf.
% Der Cache spart die zweite Zustandssimulation.
    persistent zLast Xlast
    if ~isempty(zLast) && isequal(zLast, z)
        X = Xlast; return;
    end
    u = build_u(z, VP);
    X = sim_core(VP.t_con, VP.x0(:), u, VP.p, VP.DOTstern, VP.Probe, ...
                 @Modell3_woEtOH_10p, VP.nx, 0);
    zLast = z;  Xlast = X;
end


function [FM, X, t_sim] = fim_design(z, VP)
% Simuliert das geplante Experiment mit Sensitivitaeten und bildet die FIM.
    nx = VP.nx;   np = VP.np;
    u     = build_u(z, VP);
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
            dhdx    = Modell3_dmgldx(xk);
            Sy(k,:) = dhdx(iOut,:) * XP_k;
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
            % werden zwei Matrizen unterschiedlicher Konvention addiert.
            keep = subsample_idx(kan{i,2}, VP.dt_base);
            if nnz(keep) < 3, continue; end
            dS = diff(Sy(keep,:), 1, 1);
            vk = v(keep);
            vd = vk(2:end) + vk(1:end-1);
            for k = 1:size(dS,1)
                FM = FM + (dS(k,:).' * dS(k,:)) / vd(k);
            end
        else
            for k = 1:n
                FM = FM + (Sy(k,:).' * Sy(k,:)) / v(k);
            end
        end
    end
end


function [uG, uA, uP] = unpack_u(z, VP)
% Normierte Variablen -> physikalische Feedraten auf dem feinen Raster.
    nG = VP.nu;  nA = VP.nuA;
    zG = z(1:nG);
    zA = z(nG+1 : nG+nA);
    zP = z(nG+nA+1 : end);
    uG = VP.uGlc_max * zG(:).';
    uA = VP.uAm_max  * reshape(zA(VP.mapA), 1, []);   % grob -> fein
    uP = VP.uPh_max  * reshape(zP(VP.mapP), 1, []);
end


function u = build_u(z, VP)
% Baut die u-Matrix. Letzter Wert wird gehalten (Nullter-Ordnung-Halt).
    [uG, uA, uP] = unpack_u(z, VP);
    u = zeros(10, numel(VP.tu));
    u(1,:) = VP.tu;
    u(2,:) = [uA uA(end)];   u(3,:) = VP.cAm_in;
    u(4,:) = [uP uP(end)];   u(5,:) = VP.cPh_in;
    u(6,:) = [uG uG(end)];   u(7,:) = VP.cGlc_in;
end


function X = sim_core(t, x0, u, p, DOTstern, Probe, rhs, nx, np)
% Segmentweise Integration mit Probenahme-Spruengen. Segmentgrenzen nur an
% echten Unstetigkeiten (Feedwechsel, Probenahme) -- sonst haengt das
% Ergebnis vom Auswerteraster ab.
%
% NonNegative nur fuer den reinen Zustandslauf. Im Sensitivitaetslauf waere
% die Klammerung nicht differenzierbar; dort haelt die Nebenbedingung
% cAm/cPh >= cmin die Zustaende im gueltigen Bereich.
    t  = t(:);
    tp = Probe.BatchAge(:);  vp = Probe.Volumen(:);
    tf = u(1,:).';           t0 = u(1,1);
    I  = tp >= t0 & tp <= max(t);  tp = tp(I);  vp = vp(I);
    tf = tf(tf >= t0 & tf <= max(t));
    tb = unique([t0; tf; tp; max(t)]);

    x = x0(:);
    X = zeros(numel(t), numel(x0));
    if np == 0
        opt = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.5, ...
                     'NonNegative',1:nx);
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
        for i = find(sel).'
            X(i,:) = Xs(find(tout == t(i), 1), :);
        end
        x = Xs(end,:).';
        for j = find(tp == te).', x = probe_any(x, vp(j), nx, np); end
    end
end


function xe = probe_any(xe, Vp, nx, np)
% Probenahme auf Zustand und -- falls vorhanden -- Sensitivitaeten:
% x+ = g(x), also XP+ = (dg/dx) * XP.
    x = xe(1:nx);
    if np == 0
        xe = probe_m3(x, Vp);
        return
    end
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


function [relStd, CN, CV, CVabs, Fn] = analyse_fim(FM, p, iFree)
% Relative Standardabweichungen und Kondition der NORMIERten FIM. Durch die
% Normierung sind die Diagonalelemente von CV bereits relative Varianzen.
    pf = p(iFree);
    Fn = diag(pf) * FM(iFree,iFree) * diag(pf);
    if rcond(Fn) < 1e-14
        CV = pinv(Fn);
    else
        CV = Fn \ eye(numel(iFree));
    end
    relStd = sqrt(abs(diag(CV))).';
    CVabs  = diag(pf) * CV * diag(pf);
    en = sort(eig((Fn+Fn.')/2),'descend');
    CN = en(1)/max(en(end), eps);
end


function print_worst_dir(Fn, nam)
% Eigenvektor zum kleinsten Eigenwert: welche PARAMETERKOMBINATION ist am
% schlechtesten bestimmt? Aussagekraeftiger als jede skalare Kennzahl.
    [EV, EW] = eig((Fn+Fn.')/2);
    [~, imin] = min(diag(EW));
    v = EV(:,imin);
    [~, ord] = sort(abs(v), 'descend');
    fprintf('Schlechteste Richtung: ');
    for k = 1:min(3, numel(v))
        fprintf('%+.2f*%s ', v(ord(k)), nam{ord(k)});
    end
    fprintf('\n');
end


function keep = subsample_idx(t, dt_min)
% Waehlt Punkte mit mindestens dt_min Abstand (erster Punkt immer dabei).
    keep = false(numel(t),1);  last = -inf;
    for i = 1:numel(t)
        if t(i) - last >= dt_min, keep(i) = true;  last = t(i); end
    end
end
