% main_04_parameter_fitting_LHS.m
%
% Parameteridentifikation Modell 1 (Batch) und Modell 2 (Batch + Ammonium
% + Base). Gefittet wird auf RamScDef03, validiert auf RamScDef04.
%
% Kernentscheidungen (Begruendung jeweils am Block):
%   - Mittelung pro Messkanal statt Aufsummieren
%   - Base wird in Modell 2 als Inkrement bewertet
%   - KLa ist fixiert -> gefittet wird YXO_eff = KLa*YXO
%   - Multistart per LHS

clear; clc; close all;

projectRoot = pwd;
load(fullfile(projectRoot,'Daten','Daten_Processed','Processed_Batch_Data.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
rehash; clear Modell1 Modell2

kinetic = 3;      % 3 = Monod
nCutDOT = 1;      % erste DOT-Punkte nicht fitten (Einschwingen des Sensors)

N_lhs = 100;      % billige Screening-Punkte
K_opt = 3;        % davon die besten -> teurer fmincon-Start
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% Messkanal -> Zustandsindex. Die Reihenfolge ist zugleich die
% Zustandsreihenfolge, daraus wird x0 gebaut.
spec1 = {'Biomasse',1; 'Glucose',2; 'O2',3};
spec2 = {'Biomasse',1; 'Glucose',2; 'Ammonium',3; 'Base',4; 'O2',5};

% Rechte Seiten als Handles, damit beide Modelle dieselben Hilfsfunktionen
% benutzen koennen.
RHS1 = @(t,x,p,DS) Modell1(t,x,p,kinetic,true,DS);
RHS2 = @(t,x,p,DS) Modell2(t,x,p,kinetic,DS);


%% ======================================================================
%% MODELL 1
%% ======================================================================
% Parameter: [mumax, KS, YXS, YXO_eff, KLa]
%
% KLa wird NICHT gefittet, sondern auf 1 gesetzt. DOT ist quasistationaer
% (tau = 1/KLa im Sekundenbereich, Abtastung in Minuten), deshalb bestimmt
% das DOT-Signal nur das PRODUKT KLa*YXO. KLa := 1 ist also keine Annahme,
% sondern eine Reparametrisierung: gefittet wird YXO_eff = KLa*YXO.
% pLB == pUB haelt den Wert in fmincon fest.
KLa_fix = 1.0;
p0  = [0.3;  0.5;  0.15; 300;  KLa_fix];
pLB = [0.01; 0.01; 0.01; 1;    KLa_fix];
pUB = [1.0;  500;  1.0;  1e5;  KLa_fix];

[TrainFit1, x0_tr1] = prep(TrainData, spec1, nCutDOT);
obj1 = @(p) wls_error(p, x0_tr1, TrainFit1, RHS1, spec1);
[p_opt, fval] = lhs_multistart(obj1, p0, pLB, pUB, options, N_lhs, K_opt);

fprintf('\n--- Modell1: Training RamScDef03 ---\n');
fprintf('WLS-Fehler (Training): %.4f\n', fval);
fprintf('mu_max   = %.4f 1/h\n', p_opt(1));
fprintf('K_S      = %.4f g/L\n', p_opt(2));
fprintf('Y_XS     = %.4f g/g\n', p_opt(3));
fprintf('Y_XO_eff = %.4f   (= KLa*YXO, KLa fixiert auf %.1f)\n', p_opt(4), p_opt(5));

[ValFit1, x0_val1] = prep(ValData, spec1, nCutDOT);
fprintf('WLS-Fehler (Validierung RamScDef04): %.4f\n', ...
        wls_error(p_opt, x0_val1, ValFit1, RHS1, spec1));

save(fullfile(projectRoot,'Daten','p_opt','p_opt.mat'), 'p_opt');


%% ======================================================================
%% MODELL 2
%% ======================================================================
% Parameter: [mumax, KS, YXS, Y_Bam, Y_AmX, YXO_eff, KLa]
p0_m2  = [0.3;  0.5;  0.15; 1.0;  0.05;  300;  KLa_fix];
pLB_m2 = [0.01; 0.01; 0.01; 0.01; 0.001; 1;    KLa_fix];
pUB_m2 = [1.0;  500;  1.0;  10;   1.0;   1e5;  KLa_fix];

[TrainFit2, x0_tr2] = prep(TrainData, spec2, nCutDOT);
obj2 = @(p) wls_error(p, x0_tr2, TrainFit2, RHS2, spec2);
[p_opt_m2, fval_m2] = lhs_multistart(obj2, p0_m2, pLB_m2, pUB_m2, options, N_lhs, K_opt);

fprintf('\n--- Modell2: Training RamScDef03 ---\n');
fprintf('WLS-Fehler (Training): %.4f\n', fval_m2);
fprintf('mu_max   = %.4f 1/h\n', p_opt_m2(1));
fprintf('K_S      = %.4f g/L\n', p_opt_m2(2));
fprintf('Y_XS     = %.4f g/g\n', p_opt_m2(3));
fprintf('Y_Bam    = %.4f\n',     p_opt_m2(4));
fprintf('Y_AmX    = %.4f\n',     p_opt_m2(5));
fprintf('Y_XO_eff = %.4f   (= KLa*YXO, KLa fixiert auf %.1f)\n', p_opt_m2(6), p_opt_m2(7));

[ValFit2, x0_val2] = prep(ValData, spec2, nCutDOT);
fprintf('WLS-Fehler (Validierung RamScDef04): %.4f\n', ...
        wls_error(p_opt_m2, x0_val2, ValFit2, RHS2, spec2));

save(fullfile(projectRoot,'Daten','p_opt','p_opt_Modell2.mat'), 'p_opt_m2');


%% ======================================================================
%% Beitrag pro Kanal
%% ======================================================================
% Zeigt, welcher Kanal das Guetefunktional treibt. J/nch ist der mittlere
% quadrierte Residuenwert pro Kanal -- bei korrektem Modell und korrekter
% Varianz waere er ~1.
[J,c,n] = wls_error(p_opt,    x0_tr1,  TrainFit1, RHS1, spec1);
fprintf('\n-- Modell1 Training --\n');   print_channels(spec1(:,1), c, J, n);
[J,c,n] = wls_error(p_opt,    x0_val1, ValFit1,   RHS1, spec1);
fprintf('-- Modell1 Validierung --\n');  print_channels(spec1(:,1), c, J, n);
[J,c,n] = wls_error(p_opt_m2, x0_tr2,  TrainFit2, RHS2, spec2);
fprintf('\n-- Modell2 Training --\n');   print_channels(spec2(:,1), c, J, n);
[J,c,n] = wls_error(p_opt_m2, x0_val2, ValFit2,   RHS2, spec2);
fprintf('-- Modell2 Validierung --\n');  print_channels(spec2(:,1), c, J, n);


%% ======================================================================
%% Plots und Abbildungen speichern
%% ======================================================================
ylab1 = {'c_X (g/L)','c_{Glc} (g/L)','DOT (%)'};
ylab2 = {'c_X (g/L)','c_{Glc} (g/L)','c_{Am} (g/L)','m_B (L)','DOT (%)'};

plot_experiment(TrainData, p_opt,    RHS1, spec1, ylab1, nCutDOT, 'Modell1 | Training RamScDef03');
plot_experiment(ValData,   p_opt,    RHS1, spec1, ylab1, nCutDOT, 'Modell1 | Validierung RamScDef04');
plot_experiment(TrainData, p_opt_m2, RHS2, spec2, ylab2, nCutDOT, 'Modell2 | Training RamScDef03');
plot_experiment(ValData,   p_opt_m2, RHS2, spec2, ylab2, nCutDOT, 'Modell2 | Validierung RamScDef04');

bildordner = fullfile(projectRoot,'Bilder','Batchmodelle');
figs = findobj('Type','figure');
[~, ord] = sort([figs.Number]);
for k = ord(:).'
    save_fig(figs(k), figs(k).Name, bildordner);
end


%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function [D, x0] = prep(Data, spec, nCutDOT)
% Schneidet die ersten DOT-Punkte weg und baut x0 aus den ersten
% Messwerten. x0(DOT) muss mitwandern, sonst startet die Simulation
% woanders als die verbleibenden Daten.
    D = Data;
    if nCutDOT >= 1 && numel(D.O2.t) > nCutDOT
        k = (nCutDOT+1):numel(D.O2.t);
        D.O2.t = D.O2.t(k);  D.O2.y = D.O2.y(k);  D.O2.var = D.O2.var(k);
    end
    x0 = cellfun(@(f) D.(f).y(1), spec(:,1));
end


function [p_opt, fval] = lhs_multistart(obj_fun, p0, pLB, pUB, options, N_lhs, K_opt)
% Billiges Screening ueber viele LHS-Punkte, dann fmincon von p0 und den
% K_opt besten. Log-Skalierung, weil die Parameter mehrere Groessen-
% ordnungen umspannen -- linear laegen fast alle Punkte im oberen Bereich.
    p0r = p0(:).';  pLBr = pLB(:).';  pUBr = pUB(:).';

    L = lhsdesign(N_lhs, numel(p0r));
    P = 10.^(log10(pLBr) + L .* (log10(pUBr) - log10(pLBr)));
    fixed = (pLBr == pUBr);                       % z.B. KLa
    P(:, fixed) = repmat(pLBr(fixed), N_lhs, 1);

    fprintf('LHS-Screening ueber %d Punkte ...\n', N_lhs);
    Jscreen = inf(N_lhs,1);
    for k = 1:N_lhs
        try, Jscreen(k) = obj_fun(P(k,:)); catch, end
    end
    [~, order] = sort(Jscreen);
    Pstart = [p0r; P(order(1:min(K_opt,N_lhs)), :)];

    p_opt = p0r;  fval = inf;
    for k = 1:size(Pstart,1)
        fprintf('--- fmincon Start %d/%d ---\n', k, size(Pstart,1));
        try
            [pk, Jk] = fmincon(obj_fun, Pstart(k,:), [], [], [], [], ...
                               pLBr, pUBr, [], options);
            if Jk < fval, fval = Jk;  p_opt = pk; end
        catch ME
            fprintf('  fehlgeschlagen: %s\n', ME.message);
        end
    end
    if ~isfinite(fval)
        error('lhs_multistart: kein Start erfolgreich.');
    end
    if iscolumn(p0), p_opt = p_opt(:); end
end


function [J, contrib, nch] = wls_error(p, x0, D, RHS, spec)
% Guetefunktional fuer ein Experiment.
%   Residuum:   r = (y - y_sim)/sigma
%   Gewichtung: Mittelwert pro Kanal. Ohne das wuerde DOT mit hunderten
%               Punkten das Funktional dominieren, waehrend Biomasse nur
%               ~8 Punkte hat. Dadurch ist J/nch der mittlere quadrierte
%               Residuenwert pro Kanal.
%
%   Base ist ein KUMULATIVES Integral: benachbarte Punkte teilen ihre
%   Integrationsgeschichte, ihre Fehler sind also nicht unabhaengig, und
%   sigma beschreibt nur das Rauschen einer einzelnen Zugabe. Bewertet
%   werden deshalb die Inkremente y_k - y_{k-1} -- genau die Groesse, die
%   Y_Bam parametrisiert.

    nK       = size(spec,1);
    contrib  = nan(1,nK);
    nch      = nK;
    DOTstern = max(D.O2.y);
    dt_min   = 1.0;    % h, Ausduennung der Base vor dem Differenzieren

    % Ein ODE-Lauf fuer alle Kanaele, ausgewertet an allen Messzeiten
    t_all = [];
    for i = 1:nK, t_all = [t_all; D.(spec{i,1}).t(:)]; end %#ok<AGROW>
    t_all = unique(t_all);

    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    try
        [~, X] = ode15s(@(t,x) RHS(t,x,p,DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:))), J = 1e8; return; end

    J = 0;
    for i = 1:nK
        name = spec{i,1};
        mess = D.(name);
        [~, iT] = ismember(mess.t(:), t_all);
        y_sim = X(iT, spec{i,2});

        if strcmp(name,'Base')
            % Bei ~0.2 h Abstand waere das Inkrement kleiner als sein
            % eigenes Rauschen -> vorher ausduennen.
            keep = subsample_idx(mess.t, dt_min);
            if nnz(keep) < 3, contrib(i) = 0; continue; end
            yB = mess.y(keep);  vB = max(mess.var(keep), eps);
            r  = (diff(yB) - diff(y_sim(keep))) ./ sqrt(vB(2:end) + vB(1:end-1));
        else
            r = (mess.y(:) - y_sim) ./ sqrt(max(mess.var(:), eps));
        end

        contrib(i) = mean(r.^2);
        J          = J + contrib(i);
    end
    if ~isfinite(J), J = 1e8; end
end


function keep = subsample_idx(t, dt_min)
% Waehlt Punkte mit mindestens dt_min Abstand (erster Punkt immer dabei).
    keep = false(numel(t),1);  last = -inf;
    for i = 1:numel(t)
        if t(i) - last >= dt_min, keep(i) = true;  last = t(i); end
    end
end


function print_channels(namen, contrib, J, nch)
% Beitrag jedes Kanals zu J, absolut und in Prozent.
    fprintf('  J = %.4f | J/nch = %.4f  (Erwartung ~1)\n', J, J/max(nch,1));
    for i = 1:numel(namen)
        fprintf('    %-9s %8.4f  (%5.1f %%)\n', namen{i}, contrib(i), ...
                100*contrib(i)/max(J,eps));
    end
end


function plot_experiment(Data, p, RHS, spec, ylabs, nCutDOT, titel)
% Simulation auf feinem Zeitraster, Messpunkte mit Fehlerbalken darueber.
    [D, x0]  = prep(Data, spec, nCutDOT);
    DOTstern = max(D.O2.y);

    t_end = max(cellfun(@(f) max(Data.(f).t), spec(:,1))) + 1;
    t_sim = linspace(0, t_end, 300);
    [~, X] = ode15s(@(t,x) RHS(t,x,p,DOTstern), t_sim, x0, ...
                    odeset('RelTol',1e-6,'AbsTol',1e-8));

    nK = size(spec,1);
    figure('Name', titel, 'Position', [200 60 900 180*nK]);
    for k = 1:nK
        mess = Data.(spec{k,1});
        subplot(nK,1,k);
        errorbar(mess.t, mess.y, sqrt(mess.var), 'o', ...
                 'MarkerFaceColor','b','MarkerSize',4); hold on;
        plot(t_sim, X(:,spec{k,2}), 'LineWidth', 2);
        title(spec{k,1}); ylabel(ylabs{k});
        legend('Messung \pm \sigma','Simulation','Location','best'); grid on;
    end
    xlabel('BatchAge (h)');
    sgtitle(titel);
end


function save_fig(fig, name, ordner)
% Speichert eine Figure hell (fuer Folien) als PNG und PDF.
    if ~exist(ordner,'dir'), mkdir(ordner); end
    name = regexprep(strtrim(name), '[^\w\-]', '_');
    if isempty(name), name = sprintf('Figure_%d', fig.Number); end
    ziel = fullfile(ordner, [name '.png']);

    try                                   % ab R2025a: helles Theme
        set(fig, 'Theme', 'light');  drawnow;
    catch                                 % sonst von Hand umfaerben
        set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
        for a = findall(fig,'Type','axes').'
            set(a, 'Color','w', 'XColor','k', 'YColor','k', ...
                   'GridColor',[0.15 0.15 0.15], 'GridAlpha',0.15);
            set([a.Title a.XLabel a.YLabel], 'Color', 'k');
        end
        set(findall(fig,'Type','legend'), 'TextColor','k', 'Color','w');
        set(findall(fig,'Type','text'), 'Color', 'k');
        drawnow;
    end

    exportgraphics(fig, ziel, 'Resolution', 300);
    exportgraphics(fig, strrep(ziel,'.png','.pdf'), 'ContentType','vector');
    fprintf('  gespeichert: %s\n', ziel);
end