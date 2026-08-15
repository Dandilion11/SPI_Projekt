% main_05_uncertainty_FIM_m1.m
%
% Parameterunsicherheit Modell 1 ueber die Fisher-Informationsmatrix.
%
% Kette:  df/dx, df/dp -> Sensitivitaetsgleichungen XP -> dy/dtheta
%         -> FIM -> Kovarianz (Cramer-Rao)
%
% Damit die FIM zum Fit passt:
%   - TrainData (RamScDef03), darauf wurde gefittet
%   - nCutDOT wie im Fit: erster DOT-Punkt raus, x0(3) ist der zweite
%   - Mittelung pro Kanal (/n) wie im Guetefunktional. Ohne das zaehlt DOT
%     mit seinen 60 Punkten voll, im Fit aber nur als EIN Kanal.
%   - KLa ist fixiert, die FIM wird auf die freien Parameter eingeschraenkt
%     (ueber alle fuenf waere sie singulaer)

clear; clc; close all;

projectRoot = fileparts(mfilename('fullpath'));
projectRoot = fullfile(projectRoot,'..', '..');
DATEN = fullfile(projectRoot,'Daten');        % zentrale Pfadwurzel
addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');
rehash;

load(fullfile(DATEN,'Daten_Processed','Processed_Batch_Data.mat'));
p_opt = load(fullfile(DATEN,'p_opt','p_opt.mat')).p_opt(:);

kinetic    = 3;      % Monod
withOxygen = true;
nCutDOT    = 1;      % wie im Fit

nx = 3;  np = 5;
namen = {'mumax','KS','YXS','YXO_eff','KLa'};
iFree = [1 2 3 4];                     % KLa ist fixiert
nf    = numel(iFree);

%% 1. Daten (Trainingssatz, DOT-Cut wie im Fit) --------------------------
D = TrainData;
if numel(D.O2.t) > nCutDOT
    k = (nCutDOT+1):numel(D.O2.t);
    D.O2.t = D.O2.t(k);  D.O2.y = D.O2.y(k);  D.O2.var = D.O2.var(k);
end

x0       = [D.Biomasse.y(1); D.Glucose.y(1); D.O2.y(1)];
DOTstern = max(D.O2.y);
t_all    = unique([D.Biomasse.t(:); D.Glucose.t(:); D.O2.t(:)]);

fprintf('Frei:    %s\n', strjoin(namen(iFree), ', '));
fprintf('Fixiert: %s\n', strjoin(namen(setdiff(1:np,iFree)), ', '));

%% 2. Zustaende + Sensitivitaeten simulieren -----------------------------
x0_ext = [x0; zeros(nx*np,1)];                 % XP(0) = 0
opts   = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~, X_ext] = ode15s(@(t,x) Modell1_XP_Monod(t, nx, np, x, p_opt, ...
                    withOxygen, DOTstern), t_all, x0_ext, opts);

if size(X_ext,1) ~= numel(t_all)
    error('Integration unvollstaendig -- FIM nicht berechenbar.');
end

%% 3. FIM aufsummieren ---------------------------------------------------
% {Zeitvektor, Varianz, Zustandsindex, Name}
kanal = { D.Biomasse.t, D.Biomasse.var, 1, 'Biomasse'; ...
          D.Glucose.t,  D.Glucose.var,  2, 'Glucose';  ...
          D.O2.t,       D.O2.var,       3, 'DOT'      };

FM = zeros(np, np);
for c = 1:size(kanal,1)
    [tf, idx] = ismember(kanal{c,1}(:), t_all);
    if any(~tf), error('Messzeit fuer %s nicht gefunden.', kanal{c,4}); end
    vc = kanal{c,2}(:);
    n  = numel(idx);
    for k = 1:n
        XP_k = reshape(X_ext(idx(k), nx+1:end), nx, np);
        s    = XP_k(kanal{c,3},:).';                 % dy/dtheta  (np x 1)
        FM   = FM + (s*s.') / max(vc(k),eps); % / n;    % /n = Mittelung
    end
end

%% 4. Auf die freien Parameter einschraenken -----------------------------
FM_f = FM(iFree, iFree);
p_f  = p_opt(iFree);

% Normiert (dimensionslos) -> Konditionszahl ueber die Modelle vergleichbar
FMn = diag(p_f) * FM_f * diag(p_f);
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
% Werte nahe +-1 heissen, dass nur ihr Produkt bzw. Verhaeltnis festliegt.
Corr = diag(1./max(sd,eps)) * CV * diag(1./max(sd,eps));
fprintf('\n--- Korrelationsmatrix ---\n');
fprintf('%10s',''); fprintf('%10s', namen{iFree}); fprintf('\n');
for i = 1:nf
    fprintf('%10s', namen{iFree(i)});
    fprintf('%10.2f', Corr(i,:));  fprintf('\n');
end

%% 6. Ellipsoid ----------------------------------------------------------
idx2 = [1 2];      % mumax, KS, YXS
figure('Name','Parameterunsicherheit Modell 1 (FIM)');
plot_gaussian_ellipsoid(p_f(idx2), CV(idx2,idx2), 1);
xlabel('\mu_{max}'); ylabel('K_S');
title(sprintf('1\\sigma-Ellipsoid Modell 1 (Kondition = %.2e)', en(1)/en(end)));
grid on;

if ~exist(fullfile(DATEN,'FIM'),'dir'), mkdir(fullfile(DATEN,'FIM')); end
save(fullfile(DATEN,'FIM','FIM_Modell1.mat'), ...
     'FM','FM_f','CV','sd','relStd','Corr','iFree','p_opt');
fprintf('\nGespeichert: FIM_Modell1.mat\n');

build_heatmap(Corr, namen(iFree))

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