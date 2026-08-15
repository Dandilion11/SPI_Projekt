% main_05_uncertainty_FIM_m2.m
%
% Parameterunsicherheit Modell 2 ueber die Fisher-Informationsmatrix.
% Aufbau wie Modell 1, zusaetzlich:
%   - Ammonium und Base als weitere Kanaele
%   - Base wird als INKREMENT bewertet, wie im Guetefunktional: m_B ist ein
%     kumulatives Integral, seine Residuen sind autokorreliert. Messgroesse
%     ist y_k - y_{k-1}, Sensitivitaet s_k - s_{k-1}, Varianz
%     sigma_k^2 + sigma_{k-1}^2.

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fullfile(projectRoot,'..', '..');
DATEN = fullfile(projectRoot,'Daten');        % zentrale Pfadwurzel
addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');
rehash;

load(fullfile(DATEN,'Daten_Processed','Processed_Batch_Data.mat'));
p_opt_m2 = load(fullfile(DATEN,'p_opt','p_opt_Modell2.mat')).p_opt_m2(:);

kinetic = 3;         % Monod
nCutDOT = 1;         % wie im Fit
dt_min  = 1.0;       % h, Ausduennung der Base vor dem Differenzieren

nx = 5;  np = 7;
namen = {'mumax','KS','YXS','YBam','YAmX','YXO_eff','KLa'};
iFree = [1 2 3 4 5 6];                 % KLa ist fixiert
nf    = numel(iFree);

%% 1. Daten --------------------------------------------------------------
D = TrainData;
if numel(D.O2.t) > nCutDOT
    k = (nCutDOT+1):numel(D.O2.t);
    D.O2.t = D.O2.t(k);  D.O2.y = D.O2.y(k);  D.O2.var = D.O2.var(k);
end

x0 = [D.Biomasse.y(1); D.Glucose.y(1); D.Ammonium.y(1); ...
      D.Base.y(1);     D.O2.y(1)];
DOTstern = max(D.O2.y);
t_all = unique([D.Biomasse.t(:); D.Glucose.t(:); D.Ammonium.t(:); ...
                D.Base.t(:);     D.O2.t(:)]);

fprintf('Frei:    %s\n', strjoin(namen(iFree), ', '));
fprintf('Fixiert: %s\n', strjoin(namen(setdiff(1:np,iFree)), ', '));

%% 2. Zustaende + Sensitivitaeten simulieren -----------------------------
x0_ext = [x0; zeros(nx*np,1)];                 % XP(0) = 0
opts   = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~, X_ext] = ode15s(@(t,x) Modell2_XP_Monod(t, x, p_opt_m2, DOTstern), ...
                    t_all, x0_ext, opts);

if size(X_ext,1) ~= numel(t_all)
    error('Integration unvollstaendig -- FIM nicht berechenbar.');
end

%% 3. FIM aufsummieren ---------------------------------------------------
% {Zeitvektor, Varianz, Zustandsindex, Name}
kanal = { D.Biomasse.t, D.Biomasse.var, 1, 'Biomasse'; ...
          D.Glucose.t,  D.Glucose.var,  2, 'Glucose';  ...
          D.Ammonium.t, D.Ammonium.var, 3, 'Ammonium'; ...
          D.Base.t,     D.Base.var,     4, 'Base';     ...
          D.O2.t,       D.O2.var,       5, 'DOT'       };
wsig_m2 = [1; 1; 1; 0.15; 0.0001];

FM = zeros(np, np);

for c = 1:size(kanal,1)
    tc = kanal{c,1}(:);  vc = kanal{c,2}(:);  nmc = kanal{c,4};
    [tf, idx] = ismember(tc, t_all);
    if any(~tf), error('Messzeit fuer %s nicht gefunden.', nmc); end

    % Ausgangssensitivitaeten an allen Messpunkten dieses Kanals
    Sy = zeros(numel(idx), np);
    for k = 1:numel(idx)
        XP_k    = reshape(X_ext(idx(k), nx+1:end), nx, np);
        Sy(k,:) = XP_k(kanal{c,3},:);
    end
    FC = zeros(np, np);
    if strcmp(nmc,'Base')
        keep = subsample_idx(tc, dt_min);
        if nnz(keep) < 3
            fprintf('[HINWEIS] Base: nur %d Punkte -- uebersprungen.\n', nnz(keep));
            continue
        end
        dS = diff(Sy(keep,:), 1, 1);
        vk = max(vc(keep), eps);
        vd = vk(2:end) + vk(1:end-1);
        n  = size(dS,1);
        for k = 1:n
            FC = FC + (dS(k,:).' * dS(k,:)) / vd(k); % / n;
        end
        FM = FM + wsig_m2(c) * (FC / n);
        fprintf('Base: %d Inkremente verwendet\n', n);
    else
        n = size(Sy,1);
        for k = 1:n
            FC = FC + (Sy(k,:).' * Sy(k,:)) / max(vc(k),eps); % / n;
        end
        FM = FM + wsig_m2(c) * FC / n;
    end
end

%% 4. Auf die freien Parameter einschraenken -----------------------------
FM_f = FM(iFree, iFree);
p_f  = p_opt_m2(iFree);

FMn = diag(p_f) * FM_f * diag(p_f);       % dimensionslos, modellvergleichbar
en  = sort(eig(FMn), 'descend');

fprintf('\n--- Kondition der FIM ---\n');
fprintf('Rang  = %d von %d | rcond = %.3e\n', rank(FM_f), nf, rcond(FM_f));
fprintf('Eigenwerte (normiert): %s\n', mat2str(en(:).', 3));
fprintf('Kondition (normiert)  = %.3e\n', en(1)/en(end));

if rank(FM_f) < nf
    warning('FIM singulaer -- mindestens ein Parameter nicht bestimmbar.');
end

%% 5. Kovarianz und Unsicherheiten ---------------------------------------
if rcond(FM_f) < 1e-12
    fprintf('[HINWEIS] rcond < 1e-12 -> Pseudoinverse.\n');
    CV = pinv(FM_f);
else
    CV = FM_f \ eye(nf);
end

sd     = sqrt(abs(diag(CV)));
relStd = sd ./ abs(p_f);

fprintf('\n--- Parameterunsicherheiten (FIM, 1 sigma) ---\n');
fprintf('%-10s %12s %12s %10s\n','Param','Wert','StdAbw','rel. [%]');
for i = 1:nf
    fprintf('%-10s %12.4f %12.4g %10.1f\n', ...
            namen{iFree(i)}, p_f(i), sd(i), 100*relStd(i));
end

% Korrelation: welche Parameter bestimmt die Messung nur als Kombination?
% Bei Modell 2 gehen YBam und YAmX nur als PRODUKT in die Basebildung ein,
% Ammonium liegt am Rauschgrund -> starke Antikorrelation zu erwarten.
Corr = diag(1./max(sd,eps)) * CV * diag(1./max(sd,eps));
fprintf('\n--- Korrelationsmatrix ---\n');
fprintf('%10s',''); fprintf('%10s', namen{iFree}); fprintf('\n');
for i = 1:nf
    fprintf('%10s', namen{iFree(i)});
    fprintf('%10.2f', Corr(i,:));  fprintf('\n');
end

%% 6. Ellipsoid ----------------------------------------------------------
idx2 = [1 2];      % mumax, KS, YXS
figure('Name','Parameterunsicherheit Modell 2 (FIM)');
plot_gaussian_ellipsoid(p_f(idx2), CV(idx2,idx2), 1);
xlabel('\mu_{max}'); ylabel('K_S');
title(sprintf('1\\sigma-Ellipsoid Modell 2 (Kondition = %.2e)', en(1)/en(end)));
grid on;

if ~exist(fullfile(DATEN,'FIM'),'dir'), mkdir(fullfile(DATEN,'FIM')); end
save(fullfile(DATEN,'FIM','FIM_Modell2.mat'), ...
     'FM','FM_f','CV','sd','relStd','Corr','iFree','p_opt_m2');
fprintf('\nGespeichert: FIM_Modell2.mat\n');
build_heatmap(Corr, namen(iFree))

%% ======================================================================
function keep = subsample_idx(t, dt_min)
% Waehlt Punkte mit mindestens dt_min Abstand (erster Punkt immer dabei).
    keep = false(numel(t),1);  last = -inf;
    for i = 1:numel(t)
        if t(i) - last >= dt_min, keep(i) = true;  last = t(i); end
    end
end
function build_heatmap(corr, params)
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
end