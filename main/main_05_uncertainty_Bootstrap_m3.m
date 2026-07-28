% BISHER NUR MIT AI ERSTELLT GEGENCHECKEN !!!
% 
% main_05_uncertainty_Bootstrap_m3.m
%
% Bootstrap-Parameterunsicherheit fuer Modell3_woEtOH (Fed-Batch).
% Gleiche Struktur wie die m1/m2-Bootstraps: Referenzsimulation mit p_opt,
% Pseudo-Messdaten (Referenz + Messrauschen) erzeugen, neu fitten, K-fach
% wiederholen, dann Statistik / Histogramme / Kovarianzellipsoid.
%
% Unterschiede zu m1/m2 (Fed-Batch-Besonderheiten):
%  - Modell3 hat die u-Matrix (Feeds) und eine FESTE Anfangsbedingung x0
%    (V, Massen, DOT). x0 wird daher NICHT je Lauf neu aus Messwerten
%    gebildet, sondern konstant gehalten (Data.x0).
%  - Die Konzentrations-Messgroessen (Biomasse/Glucose/Ammonium/Phosphat)
%    sind Massen/V; Base (mB) und DOT sind direkte Zustaende.
%  - Es wird EINMAL ueber das vereinigte Zeitraster simuliert und den
%    jeweiligen Messzeitpunkten zugeordnet (wie in wls_error_m3).
close all; clear all; clc;

%% M3.1 Laden
scriptDir   = pwd;
projectRoot = fileparts(scriptDir);
load(fullfile(scriptDir,'..','Daten','Daten_Processed','Processed_FedBatch_Modell3.mat'));
tmp3  = load(fullfile(scriptDir,'..','Daten','p_opt_Modell3_woEtOH.mat'));
p_opt_m3 = tmp3.p_opt(:);   % [mumax, KS, YXS, YAmX, YPhX, YB_Am, KLa, YXO]

addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');

%% M3.2 Daten & Konfiguration
Data = TrainData;
u    = Data.u;
x0   = Data.x0;                 % feste Anfangsbedingung [V; mX; mGlc; mAm; mPh; mB; DOT]
DOTstern = max(Data.O2.y);      % wie im Fitting (wls_error_m3)

np3    = 8;
namen3 = {'mumax','KS','YXS','YAmX','YPhX','YB_Am','KLa','YXO'};
lb3    = [0.01, 0.01, 0.05, 0.001, 0.001, 0.1,  10, 0.1];
ub3    = [1.00, 5.00, 1.00, 1.000, 1.000, 5.0, 800, 5.0];

% Messgroessen: {Struct, Zustandsindex, divByV}
%   divByV = true  -> Konzentration = Masse / V
%   divByV = false -> direkter Zustand (mB, DOT)
M = { Data.Biomasse, 2, true;  ...
      Data.Glucose,  3, true;  ...
      Data.Ammonium, 4, true;  ...
      Data.Phosphat, 5, true;  ...
      Data.Base,     6, false; ...
      Data.O2,       7, false  };
nSig = size(M,1);

%% M3.3 Vereinigtes Zeitraster (inkl. Startzeit t0) + Zuordnungsindizes
t0    = u(1,1);
t_all = t0;
for i = 1:nSig
    t_all = [t_all; M{i,1}.t(:)];
end
t_all = unique(t_all);
if t_all(1) > t0, t_all = [t0; t_all]; end

% je Signal die Indizes im t_all-Raster
iT = cell(nSig,1);
for i = 1:nSig
    [~, iT{i}] = ismember(M{i,1}.t(:), t_all);
end

%% M3.4 Referenzsimulation mit p_opt  ->  Basis fuer Pseudo-Messdaten
opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~, Xref] = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p_opt_m3,DOTstern), t_all, x0, opts_ode);
Vref = Xref(:,1);

y_ref = cell(nSig,1);
for i = 1:nSig
    yi = Xref(iT{i}, M{i,2});
    if M{i,3}, yi = yi ./ Vref(iT{i}); end   % Konzentration = Masse/V
    y_ref{i} = yi;
end

%% M3.5 Bootstrap-Schleife
K3        = 50;                 % Anzahl Laeufe
p_boot_m3 = zeros(K3, np3);
rng(42);
opts_opt = optimoptions('fmincon','Display','off','Algorithm','sqp');
fprintf('\nBootstrap Modell3_woEtOH (%d Laeufe)...\n', K3);

for k = 1:K3
    % --- Pseudo-Messdaten erzeugen (Referenz + Messrauschen) ---
    y_boot = cell(nSig,1);
    for i = 1:nSig
        vi = M{i,1}.var(:);
        y_boot{i} = y_ref{i} + sqrt(vi) .* randn(size(y_ref{i}));
    end

    % --- Parameterschaetzung auf Pseudo-Daten (Start = p_opt) ---
    obj_k = @(p) boot_wls_m3(p, x0, u, t_all, iT, M, y_boot, DOTstern);
    try
        p_boot_m3(k,:) = fmincon(obj_k, p_opt_m3, [],[],[],[], lb3, ub3, [], opts_opt);
    catch
        p_boot_m3(k,:) = NaN;
    end

    if mod(k,25)==0, fprintf('  %d/%d\n',k,K3); end
end

%% M3.6 Ungueltige Laeufe entfernen
valid3    = ~any(isnan(p_boot_m3),2);
p_boot_m3 = p_boot_m3(valid3,:);
fprintf('Gueltige Laeufe: %d/%d\n', sum(valid3), K3);

%% M3.7 Statistik
p_mean_m3  = mean(p_boot_m3,1);
CV_boot_m3 = cov(p_boot_m3);

[Corr3, StdDev3, relStdDev3, ~, ~, CN3] = Parameteranalyse(CV_boot_m3, p_opt_m3(:));

fprintf('\n--- Bootstrap Parameterunsicherheiten (Modell3_woEtOH) ---\n');
fprintf('%-8s  %-10s  %-10s  %-12s  %-14s\n','Param','p_opt','p_mean','StdAbw','Rel. [%]');
for i = 1:np3
    fprintf('%-8s  %-10.4f  %-10.4f  %-12.4f  %-14.1f\n', ...
        namen3{i}, p_opt_m3(i), p_mean_m3(i), StdDev3(i), relStdDev3(i)*100);
end
fprintf('\nKorrelationsmatrix:\n'); disp(Corr3);
fprintf('Konditionszahl: %.1f\n', CN3);

%% M3.8 Histogramme
nCols3 = 4;
nRows3 = ceil(np3/nCols3);
figure('Name','Bootstrap Modell3_woEtOH – Parameterverteilungen','Position',[120 80 480*nCols3 500*nRows3]);
for i = 1:np3
    subplot(nRows3,nCols3,i);
    histogram(p_boot_m3(:,i), 25, 'Normalization','pdf'); hold on;
    xline(p_opt_m3(i),  'r-',  'LineWidth',2);
    xline(p_mean_m3(i), 'g--', 'LineWidth',1.5);
    xlabel(namen3{i}); ylabel('Dichte');
    title(sprintf('%s:  \\sigma=%.4g (%.1f%%)', namen3{i}, StdDev3(i), relStdDev3(i)*100));
    legend('Bootstrap','p_{opt}','\mu_{boot}',Location='northeast');
    grid on;
end
sgtitle(sprintf('Bootstrap-Verteilungen Modell3\\_woEtOH  (K=%d)', sum(valid3)));

%% M3.9 Kovarianzellipsoid (Wachstumsparameter [mu_max, K_S, Y_XS])
idxEll3 = [1 2 3];
figure('Name','Parameterunsicherheit Modell3_woEtOH (Bootstrap)','Position',[220 180 600 500]);
plot_gaussian_ellipsoid(p_opt_m3(idxEll3), CV_boot_m3(idxEll3,idxEll3), 1);
xlabel('\mu_{max}'); ylabel('K_S'); zlabel('Y_{XS}');
title(sprintf('1\\sigma-Ellipsoid Wachstumsparameter (Modell3\\_woEtOH, K=%d)', sum(valid3)));
grid on;


%% ---- Lokale Hilfsfunktion ----------------------------------------
function J = boot_wls_m3(p, x0, u, t_all, iT, M, y_boot, DOTstern)
    % WLS-Guete Modell3_woEtOH auf Pseudo-Daten.
    % Eine Integration ueber t_all (Start bei x0), dann Zuordnung zu den
    % Messzeitpunkten je Signal; Konzentrationen = Masse/V.
    opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
    try
        [~, X] = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p,DOTstern), t_all, x0, opts);
    catch
        J = 1e6; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:)))
        J = 1e6; return;
    end

    V = X(:,1);
    J = 0;
    for i = 1:size(M,1)
        idxState = M{i,2};  divByV = M{i,3};
        y_sim = X(iT{i}, idxState);
        if divByV, y_sim = y_sim ./ V(iT{i}); end
        vi = M{i,1}.var(:);
        J  = J + sum(((y_sim - y_boot{i}).^2) ./ vi);
    end

    if ~isfinite(J), J = 1e6; end
end
