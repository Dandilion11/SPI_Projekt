% main_05_uncertainty_Bootstrap.m
clear; clc; close all;

%% 1. Laden
scriptDir   = pwd;
projectRoot = fileparts(scriptDir);

load(fullfile(scriptDir,'..','Daten','Daten_Processed','Processed_Batch_Data.mat'));
tmp   = load(fullfile(scriptDir,'..','Daten','p_opt.mat'));
p_opt = tmp.p_opt;

addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');

%% 2. Messdaten
t_mes   = ValData.Biomasse.t;
y_bio   = ValData.Biomasse.y;
y_glc   = ValData.Glucose.y;
var_bio = ValData.Biomasse.var;
var_glc = ValData.Glucose.var;

x0         = [y_bio(1); y_glc(1)];
kinetic    = 3;
withOxygen = false;
np         = 3;
lb = [0.01, 0.01, 0.01];
ub = [1.0,  5.0,  1.0 ];

%% 3. Referenzsimulation mit p_opt  →  Basis für Pseudo-Messdaten
opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~, X_ref] = ode15s(@(t,x) Modell1(t,x,p_opt,kinetic,withOxygen), t_mes, x0, opts_ode);
y_bio_ref = X_ref(:,1);
y_glc_ref = X_ref(:,2);

%% 4. Bootstrap-Schleife
K      = 300;                   % Anzahl Läufe
p_boot = zeros(K, np);
rng(42);

opts_opt = optimoptions('fmincon','Display','off','Algorithm','sqp');
fprintf('Bootstrap (%d Läufe)...\n', K);

for k = 1:K
    % --- Pseudo-Messdaten erzeugen ---
    y_bio_k = y_bio_ref + sqrt(var_bio) .* randn(size(y_bio_ref));
    y_glc_k = y_glc_ref + sqrt(var_glc) .* randn(size(y_glc_ref));

    x0_k = max([y_bio_k(1); y_glc_k(1)], 1e-6);

    % --- Parameterschätzung auf Pseudo-Daten ---
    obj_k = @(p) boot_wls(p, x0_k, t_mes, y_bio_k, y_glc_k, ...
                           var_bio, var_glc, kinetic, withOxygen);
    try
        % Startwert = p_opt  -> schnelle Konvergenz
        p_boot(k,:) = fmincon(obj_k, p_opt, [],[],[],[], lb, ub, [], opts_opt);
    catch
        p_boot(k,:) = NaN;
    end

    if mod(k,50)==0, fprintf('  %d/%d\n',k,K); end
end

%% 5. Ungültige Läufe entfernen
valid  = ~any(isnan(p_boot),2);
p_boot = p_boot(valid,:);
fprintf('Gültige Läufe: %d/%d\n', sum(valid), K);

%% 6. Statistik
p_mean  = mean(p_boot,1);
CV_boot = cov(p_boot);

[Corr, StdDev, relStdDev, ~, ~, CN] = Parameteranalyse(CV_boot, p_opt(:));

namen = {'mumax','KS','YXS'};
fprintf('\n--- Bootstrap Parameterunsicherheiten ---\n');
fprintf('%-8s  %-10s  %-10s  %-12s  %-14s\n','Param','p_opt','p_mean','StdAbw','Rel. [%]');
for i = 1:np
    fprintf('%-8s  %-10.4f  %-10.4f  %-12.4f  %-14.1f\n', ...
        namen{i}, p_opt(i), p_mean(i), StdDev(i), relStdDev(i)*100);
end
fprintf('\nKorrelationsmatrix:\n'); disp(Corr);
fprintf('Konditionszahl: %.1f\n', CN);

%% 7. Histogramme
figure('Name','Bootstrap – Parameterverteilungen','Position',[100 100 900 350]);
for i = 1:np
    subplot(1,np,i);
    histogram(p_boot(:,i), 25, 'Normalization','pdf'); hold on;
    xline(p_opt(i),  'r-',  'LineWidth',2);
    xline(p_mean(i), 'g--', 'LineWidth',1.5);
    xlabel(namen{i}); ylabel('Dichte');
    title(sprintf('%s:  \\sigma=%.4f (%.1f%%)', namen{i}, StdDev(i), relStdDev(i)*100));
    legend('Bootstrap','p_{opt}','\mu_{boot}','Location','best');
    grid on;
end
sgtitle(sprintf('Bootstrap-Verteilungen  (K=%d)', sum(valid)));

%% 8. Kovarianzellipsoid
figure('Name','Parameterunsicherheit (Bootstrap)','Position',[200 200 600 500]);
plot_gaussian_ellipsoid(p_opt(:), CV_boot, 1);
xlabel('\mu_{max}'); ylabel('K_S'); zlabel('Y_{XS}');
title(sprintf('1\\sigma-Ellipsoid (Bootstrap, K=%d)', sum(valid)));
grid on;

%% ---- Lokale Hilfsfunktion ----------------------------------------
function J = boot_wls(p, x0, t, y_bio, y_glc, v_bio, v_glc, kin, wO2)
    try
        [~,X] = ode15s(@(t,x) Modell1(t,x,p,kin,wO2), t, x0, ...
                       odeset('RelTol',1e-4,'AbsTol',1e-6));
        J = sum(((X(:,1)-y_bio(:)).^2)./v_bio(:)) + ...
            sum(((X(:,2)-y_glc(:)).^2)./v_glc(:));
        if ~isfinite(J), J = 1e6; end
    catch
        J = 1e6;
    end
end