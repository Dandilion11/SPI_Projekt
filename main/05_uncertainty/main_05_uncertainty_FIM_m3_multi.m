% main_05_uncertainty_FIM_m3_multi.m
%
% Parameterunsicherheit Modell 3 ueber die Fisher-Informationsmatrix.
%
% Kette:  df/dx, df/dp -> Sensitivitaetsgleichungen XP -> dy/dtheta
%         -> FIM -> Kovarianz (Cramer-Rao)
%
% Aufbau wie bei Modell 1/2, mit drei Besonderheiten des Fed-Batch:
%   - Summe ueber ALLE Trainingsexperimente, so wurde auch gefittet
%   - Probenahmen: an jedem Probenzeitpunkt springt nicht nur x, sondern
%     auch XP (Kettenregel ueber die Probenahme-Abbildung)
%   - Messgroessen sind Konzentrationen, Zustaende sind Massen -> die
%     Umrechnung steckt in Modell3_dmgldx
%
% CAVEAT fuer die Interpretation: die FIM beschreibt die Streuung der
% Schaetzung unter MESSRAUSCHEN bei korrektem Modell. chi2/N_eff liegt aber
% bei ~41, der Strukturfehler (fehlende Ethanol-Kinetik) dominiert also.
% Die Intervalle sind eine untere Schranke der statistischen Praezision,
% keine Aussage ueber die Richtigkeit der Parameter.

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fullfile(projectRoot,'..', '..');
DATEN = fullfile(projectRoot,'Daten');        % zentrale Pfadwurzel
addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');
rehash;

load(fullfile(DATEN,'Daten_Processed','Processed_FedBatch_Modell3_MultiExp.mat'));
S = load(fullfile(DATEN,'p_opt','p_opt_Modell3_woEtOH_10p_multi.mat'));
p_opt = S.p_opt(:);

nx = 7;   np = 10;
namen = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO','KAm','KPh'};
iFree = [1 2 3 4 5 6 8];  % KLa, KAm, KPh waren im Fit fixiert
nf    = numel(iFree);

fprintf('Frei:    %s\n', strjoin(namen(iFree), ', '));
fprintf('Fixiert: %s\n', strjoin(namen(setdiff(1:np,iFree)), ', '));

%% 1. FIM ueber alle Trainingsexperimente aufsummieren -------------------
FM = zeros(np, np);
fprintf('\n--- Beitrag je Experiment ---\n');
for k = 1:numel(TrainSet)
    D  = TrainSet(k);
    Fk = fim_one(D, p_opt, max(D.O2.y), nx, np, S.wsig);
    FM = FM + Fk;
    fprintf('%-12s  rank = %d | trace = %.3e\n', ...
            D.name, rank(Fk(iFree,iFree)), trace(Fk(iFree,iFree)));
end

%% 2. Auf die freien Parameter einschraenken -----------------------------
FM_f = FM(iFree, iFree);
p_f  = p_opt(iFree);

FMn = diag(p_f) * FM_f * diag(p_f);       % dimensionslos, modellvergleichbar
en  = sort(eig(FMn), 'descend');

fprintf('\n--- Kondition der FIM ---\n');
fprintf('Rang  = %d von %d | rcond = %.3e\n', rank(FM_f), nf, rcond(FM_f));
fprintf('Eigenwerte (normiert): %s\n', mat2str(en(:).', 3));
fprintf('Kondition (normiert)  = %.3e\n', en(1)/en(end));

if rank(FM_f) < nf
    warning('FIM singulaer -- mindestens ein Parameter nicht bestimmbar.');
end

%% 3. Kovarianz und Unsicherheiten ---------------------------------------
if rcond(FM_f) < 1e-12
    fprintf('[HINWEIS] rcond < 1e-12 -> Pseudoinverse.\n');
    CV = pinv(FM_f);
else
    CV = FM_f \ eye(nf);
end

sd     = sqrt(abs(diag(CV)));
relStd = sd ./ abs(p_f);

fprintf('\n--- Parameterunsicherheiten (FIM, 1 sigma) ---\n');
fprintf('%-9s %12s %12s %10s\n','Param','Wert','StdAbw','rel. [%]');
for i = 1:nf
    fprintf('%-9s %12.4f %12.4g %10.1f\n', ...
            namen{iFree(i)}, p_f(i), sd(i), 100*relStd(i));
end

% Korrelation: welche Parameter bestimmt die Messung nur als Kombination?
% YXS und YAmX skalieren beide denselben Term rX*mX -> stark antikorreliert.
Corr = diag(1./max(sd,eps)) * CV * diag(1./max(sd,eps));
fprintf('\n--- Korrelationsmatrix ---\n');
fprintf('%9s',''); fprintf('%9s', namen{iFree}); fprintf('\n');
for i = 1:nf
    fprintf('%9s', namen{iFree(i)});
    fprintf('%9.2f', Corr(i,:));  fprintf('\n');
end

%% 4. Korrektur um den Strukturfehler ------------------------------------
% Die FIM unterstellt, die Residuen seien reines Messrauschen. Bei
% chi2/N_eff = c sind sie im Mittel sqrt(c) mal groesser -- die skalierte
% Spalte ist die ehrlichere Angabe fuer "wie gut kennen wir den Parameter".
chi2_N = S.fval / (nnz(S.wsig) * numel(TrainSet));
fprintf('\nchi2/N_eff = %.1f  ->  Skalierung mit sqrt = %.1f\n', ...
        chi2_N, sqrt(chi2_N));
fprintf('%-9s %12s %12s\n','Param','sd (FIM)','sd skaliert');
for i = 1:nf
    fprintf('%-9s %12.4g %12.4g\n', namen{iFree(i)}, sd(i), sd(i)*sqrt(chi2_N));
end

%% 5. Ellipsoid ----------------------------------------------------------
idx2 = [3 4];  % YXS, YAmX (Positionen innerhalb von iFree)
figure('Name','Parameterunsicherheit Modell 3 (FIM)');
plot_gaussian_ellipsoid(p_f(idx2), CV(idx2,idx2), 1);
xlabel('Y_{XS}'); ylabel('Y_{AmX}');
title(sprintf('1\\sigma-Ellipsoid Modell 3 (Kondition = %.2e)', en(1)/en(end)));
grid on;

if ~exist(fullfile(DATEN,'FIM'),'dir'), mkdir(fullfile(DATEN,'FIM')); end
save(fullfile(DATEN,'FIM','FIM_Modell3_multi.mat'), ...
     'FM','FM_f','CV','sd','relStd','Corr','iFree','p_opt','chi2_N');
fprintf('\nGespeichert: FIM_Modell3_multi.mat\n');

hFig = build_heatmap(Corr, namen(iFree));

% Heatmap speichern
bildordner = fullfile(projectRoot,'Bilder','FIM');
if ~exist(bildordner,'dir'), mkdir(bildordner); end
set(hFig, 'Color', 'w', 'InvertHardcopy', 'off');
drawnow;
try
    exportgraphics(hFig, fullfile(bildordner,'Korrelationsmatrix_M3.svg'), ...
        'ContentType', 'vector');
catch
    set(hFig, 'PaperPositionMode', 'auto');
    print(hFig, fullfile(bildordner,'Korrelationsmatrix_M3.svg'), '-dsvg');
end
fprintf('  gespeichert: Korrelationsmatrix_M3.svg\n');

%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function F = fim_one(D, p, DOTstern, nx, np, wsig)
    M  = { D.Biomasse, 1; D.Glucose, 2; D.Ammonium, 3; ...
           D.Phosphat, 4; D.Base,    5; D.O2,       6 };
    nm = {'Biomasse','Glucose','Ammonium','Phosphat','Base','DOT'};

    t_all = [];
    for i = 1:6, t_all = [t_all; M{i,1}.t(:)]; end %#ok<AGROW>
    t_all = unique([t_all; D.u(1,:).']);
    t_all = t_all(t_all >= D.u(1,1));

    X_ext = sim_xp(t_all, D.x0, D.u, p, DOTstern, D.Probe, nx, np);

    F = zeros(np, np);
    for i = 1:6
        if wsig(i) == 0, continue; end
        mess = M{i,1};
        [tf, iT] = ismember(mess.t(:), t_all);
        if any(~tf), error('Messzeitpunkt fuer %s nicht gefunden.', nm{i}); end

        Sy = zeros(numel(iT), np);
        for k = 1:numel(iT)
            xk   = X_ext(iT(k), 1:nx).';
            XP_k = reshape(X_ext(iT(k), nx+1:end), nx, np);
            dhdx = Modell3_dmgldx(xk);
            Sy(k,:) = dhdx(M{i,2},:) * XP_k;
        end

        v = max(mess.var(:), eps);
        n = size(Sy,1);
        for k = 1:n
            F = F + (wsig(i)/n) * (Sy(k,:).' * Sy(k,:)) / v(k);
        end
    end
end

function X_ext = sim_xp(t, x0, u, p, DOTstern, Probe, nx, np)
% Integriert Zustaende UND Sensitivitaeten, segmentweise zwischen den
% echten Unstetigkeiten. An einer Probenahme gilt x+ = g(x), also
% XP+ = (dg/dx) * XP.
    t  = t(:);

    tp = Probe.BatchAge(:);  vp = Probe.Volumen(:);
    tf = u(1,:).';           t0 = u(1,1);

    I  = tp >= t0 & tp <= max(t);  tp = tp(I);  vp = vp(I);
    tf = tf(tf >= t0 & tf <= max(t));
    tb = unique([t0; tf; tp; max(t)]);

    x = [x0(:); zeros(nx*np,1)];        % XP(0) = 0
    X_ext = zeros(numel(t), nx + nx*np);
    opt = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1.0);

    X_ext(t == t0,:) = repmat(x.', nnz(t == t0), 1);
    for j = find(tp == t0).', x = probe_ext(x, vp(j), nx, np); end

    for k = 2:numel(tb)
        ta = tb(k-1);  te = tb(k);
        sel  = t > ta & t <= te;
        tout = unique([ta; t(sel); te]);

        [~, Xs] = ode15s(@(tt,xx) Modell3_woEtOH_XP(tt,xx,u,p,DOTstern), ...
                         tout, x, opt);
        if numel(tout) == 2, Xs = Xs([1 end],:); end

        for i = find(sel).'
            X_ext(i,:) = Xs(find(tout == t(i), 1), :);
        end
        x = Xs(end,:).';
        for j = find(tp == te).', x = probe_ext(x, vp(j), nx, np); end
    end
end


function xe = probe_ext(xe, Vp, nx, np)
% Probenahme auf Zustand UND Sensitivitaeten anwenden. Die Jacobi-Matrix
% der Probenahme wird numerisch gebildet (probe_m3 ist kurz und glatt).
    x  = xe(1:nx);
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


function h = build_heatmap(corr, params)
n = 256;

blue = [0.0000 0.4470 0.7410];
white = [1.0000 1.0000 1.0000];
red  = [0.8500 0.3250 0.0980];

cmap1 = [linspace(blue(1),white(1),n/2)', ...
         linspace(blue(2),white(2),n/2)', ...
         linspace(blue(3),white(3),n/2)'];

cmap2 = [linspace(white(1),red(1),n/2)', ...
         linspace(white(2),red(2),n/2)', ...
         linspace(white(3),red(3),n/2)'];

cmap = [cmap1; cmap2];

figure('Color','w');

h = heatmap(params, params, corr);

h.Title = 'Korrelationsmatrix';
h.XLabel = 'Parameter';
h.YLabel = 'Parameter';

h.ColorLimits = [-1 1];
h.Colormap = cmap;
h.CellLabelFormat = '%.2f';
h.FontSize = 11;
h = gcf;
end