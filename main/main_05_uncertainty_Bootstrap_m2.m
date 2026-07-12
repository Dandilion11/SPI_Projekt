% main_05_uncertainty_Bootstrap.m

% Die DOT-Messung hat einen eigenen, viel feineren Zeitvektor
% als Biomasse/Glucose. Daher wird der O2-Anteil separat auf dem O2-Raster
% (t_o2) simuliert, Biomasse/Glucose auf ihrem gemeinsamen Raster (t_mes)
%% MODELL 2 Bootstrap
close all; clear all; clc;
scriptDir   = pwd;
projectRoot = fileparts(scriptDir);
load(fullfile(scriptDir,'..','Daten','Daten_Processed','Processed_Batch_Data.mat'));
addpath(fullfile(projectRoot,'Modelle'),'-begin');
addpath(fullfile(projectRoot,'utils'),  '-begin');
% M2.1 Modell2-Parameter laden
tmp2     = load(fullfile(scriptDir,'..','Daten','p_opt_Modell2.mat'));
p_opt_m2 = tmp2.p_opt_m2;
p_opt_m2 = p_opt_m2(:);   % als Spaltenvektor sicherstellen
 
% M2.2 Messdaten & Konfiguration
% Zwei Zeitgitter:
%   t_off (Offline): Biomasse, Glucose, Ammonium -> Zustaende 1,2,3
%   t_on  (Online) : Base, O2                    -> Zustaende 4,5
t_off = ValData.Biomasse.t;
t_on  = ValData.O2.t;
 
y_bio2 = ValData.Biomasse.y;  var_bio2 = ValData.Biomasse.var;
y_glc2 = ValData.Glucose.y;   var_glc2 = ValData.Glucose.var;
y_am   = ValData.Ammonium.y;  var_am   = ValData.Ammonium.var;
y_ba   = ValData.Base.y;      var_ba   = ValData.Base.var;
y_o22  = ValData.O2.y;        var_o22  = ValData.O2.var;
 
np2    = 7;
lb2    = [0.01, 0.01, 0.01, 0.01, 0.001, 0.01,   1];
ub2    = [1.0,  5.0,  1.0,  10,   1.0,   5.0, 1000];
namen2 = {'mumax','KS','YXS','YBam','YAmX','YXO','KLa'};
 
x0_m2 = [y_bio2(1); y_glc2(1); y_am(1); y_ba(1); y_o22(1)];
 
% M2.3 Referenzsimulation mit p_opt_m2 (je Zeitgitter einmal)
kinetic    = 3;
opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~, Xoff] = ode15s(@(t,x) Modell2(t,x,p_opt_m2,kinetic), t_off, x0_m2, opts_ode);
[~, Xon ] = ode15s(@(t,x) Modell2(t,x,p_opt_m2,kinetic), t_on , x0_m2, opts_ode);
y_bio_ref2 = Xoff(:,1);
y_glc_ref2 = Xoff(:,2);
y_am_ref   = Xoff(:,3);
y_ba_ref   = Xon(:,4);
y_o2_ref2  = Xon(:,5);
 
% M2.4 Bootstrap-Schleife
K2        = 300;
p_boot_m2 = zeros(K2, np2);
rng(42);
opts_opt = optimoptions('fmincon','Display','off','Algorithm','sqp');
fprintf('\nBootstrap Modell2 (%d Laeufe)...\n', K2);
 
for k = 1:K2
    % --- Pseudo-Messdaten erzeugen (Referenz + Messrauschen) ---
    yb  = y_bio_ref2 + sqrt(var_bio2(:)) .* randn(size(y_bio_ref2));
    yg  = y_glc_ref2 + sqrt(var_glc2(:)) .* randn(size(y_glc_ref2));
    ya  = y_am_ref   + sqrt(var_am(:))   .* randn(size(y_am_ref));
    yba = y_ba_ref   + sqrt(var_ba(:))   .* randn(size(y_ba_ref));
    yo  = y_o2_ref2  + sqrt(var_o22(:))  .* randn(size(y_o2_ref2));
 
    x0k = max([yb(1); yg(1); ya(1); yba(1); yo(1)], 1e-6);
 
    % --- Parameterschaetzung auf Pseudo-Daten (Start = p_opt_m2) ---
    obj_k = @(p) boot_wls_m2(p, x0k, t_off, t_on, yb, yg, ya, yba, yo, var_bio2, var_glc2, var_am, var_ba, var_o22, kinetic);
    try
        p_boot_m2(k,:) = fmincon(obj_k, p_opt_m2, [],[],[],[], lb2, ub2, [], opts_opt);
    catch
        p_boot_m2(k,:) = NaN;
    end
 
    if mod(k,50)==0, fprintf('  %d/%d\n',k,K2); end
end
 
% M2.5 Ungueltige Laeufe entfernen
valid2    = ~any(isnan(p_boot_m2),2);
p_boot_m2 = p_boot_m2(valid2,:);
fprintf('Gueltige Laeufe: %d/%d\n', sum(valid2), K2);
 
% M2.6 Statistik
p_mean_m2  = mean(p_boot_m2,1);
CV_boot_m2 = cov(p_boot_m2);
 
[Corr2, StdDev2, relStdDev2, ~, ~, CN2] = Parameteranalyse(CV_boot_m2, p_opt_m2(:));
 
fprintf('\n--- Bootstrap Parameterunsicherheiten (Modell2) ---\n');
fprintf('%-8s  %-10s  %-10s  %-12s  %-14s\n','Param','p_opt','p_mean','StdAbw','Rel. [%]');
for i = 1:np2
    fprintf('%-8s  %-10.4f  %-10.4f  %-12.4f  %-14.1f\n', ...
        namen2{i}, p_opt_m2(i), p_mean_m2(i), StdDev2(i), relStdDev2(i)*100);
end
fprintf('\nKorrelationsmatrix:\n'); disp(Corr2);
fprintf('Konditionszahl: %.1f\n', CN2);
 
%% M2.7 Histogramme
nCols2 = 4;
nRows2 = ceil(np2/nCols2);
figure('Name','Bootstrap Modell2 – Parameterverteilungen','Position',[120 80 500*nCols2 520*nRows2]);
for i = 1:np2
    subplot(nRows2,nCols2,i);
    histogram(p_boot_m2(:,i), 25, 'Normalization','pdf'); hold on;
    xline(p_opt_m2(i),  'r-',  'LineWidth',2);
    xline(p_mean_m2(i), 'g--', 'LineWidth',1.5);
    xlabel(namen2{i}); ylabel('Dichte');
    title(sprintf('%s:  \\sigma=%.4g (%.1f%%)', namen2{i}, StdDev2(i), relStdDev2(i)*100));
    legend('Bootstrap','p_{opt}','\mu_{boot}',Location='northeast');
    grid on;
end
sgtitle(sprintf('Bootstrap-Verteilungen Modell2  (K=%d)', sum(valid2)));
 
%% M2.8 Kovarianzellipsoid (Wachstumsparameter [mu_max, K_S, Y_XS])
idxEll2 = [1 2 3];
figure('Name','Parameterunsicherheit Modell2 (Bootstrap)','Position',[220 180 600 500]);
plot_gaussian_ellipsoid(p_opt_m2(idxEll2), CV_boot_m2(idxEll2,idxEll2), 1);
xlabel('\mu_{max}'); ylabel('K_S'); zlabel('Y_{XS}');
title(sprintf('1\\sigma-Ellipsoid Wachstumsparameter (Modell2, K=%d)', sum(valid2)));
grid on;
 

%%
function J = boot_wls_m2(p, x0, t_off, t_on, y_bio, y_glc, y_am, y_ba, y_o2, v_bio, v_glc, v_am, v_ba, v_o2, kin)
    % WLS-Guete Modell2 auf Pseudo-Daten.
    % Zwei Integrationen: Offline-Groessen (Biomasse/Glucose/Ammonium) auf
    % t_off, Online-Groessen (Base/O2) auf t_on
    opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
    try
        [~,Xo] = ode15s(@(t,x) Modell2(t,x,p,kin), t_off, x0, opts);
        [~,Xn] = ode15s(@(t,x) Modell2(t,x,p,kin), t_on , x0, opts);
        if size(Xo,1) ~= numel(t_off) || size(Xn,1) ~= numel(t_on)
            J = 1e6; return;
        end
        J = sum(((Xo(:,1)-y_bio(:)).^2)./v_bio(:)) + ...
            sum(((Xo(:,2)-y_glc(:)).^2)./v_glc(:)) + ...
            sum(((Xo(:,3)-y_am(:) ).^2)./v_am(:))  + ...
            sum(((Xn(:,4)-y_ba(:) ).^2)./v_ba(:))  + ...
            sum(((Xn(:,5)-y_o2(:) ).^2)./v_o2(:));
        if ~isfinite(J), J = 1e6; end
    catch
        J = 1e6;
    end
end



% function J = boot_wls(p, x0, t, y_bio, y_glc, v_bio, v_glc, kin, wO2)
%     try
%         [~,X] = ode15s(@(t,x) Modell1(t,x,p,kin,wO2), t, x0, ...
%                        odeset('RelTol',1e-4,'AbsTol',1e-6));
%         J = sum(((X(:,1)-y_bio(:)).^2)./v_bio(:)) + ...
%             sum(((X(:,2)-y_glc(:)).^2)./v_glc(:));
%         if ~isfinite(J), J = 1e6; end
%     catch
%         J = 1e6;
%     end
% end