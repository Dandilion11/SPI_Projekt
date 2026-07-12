% main_04_parameter_fitting.m
clear; clc; close all;

% 1. Load Preprocessed Data
projectRoot = pwd;
load(fullfile(projectRoot,'..','Daten/Daten_Processed/Processed_FedBatch_Modell3.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
rehash; clear Modell3

Data = TrainData;
u    = Data.u;
x0   = Data.x0;                % [V; mX; mGlc; mAm; mPh; mB; DOT]


% 1. Anfangswerte und Parameter
%Ausgangslage
namen = {'mumax','KS','YXS', 'YAmX','YPhX','YB_Am','KLa','YXO'};
p0 =    [0.30,   0.50, 0.15,  0.05,  0.02,  1.0,     200,  1.0];
pLB =   [0.01,   0.01, 0.05, 0.001, 0.001,  0.1,     10,   0.1];
pUB =   [1.00,   5.00, 1.00, 1.000, 1.000,  5.0,     800,  5.0];


% 2. Parameteridentifikation (WLS)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', 5000, ...
                       'FiniteDifferenceStepSize', 1e-6, ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance', 1e-10);

obj_fun = @(p) wls_error_m3(p,x0, u, Data);

[p_opt, fval] = fmincon(obj_fun, p0, [], [], [], [], ...
                        pLB, pUB, [], options);

% 3. Ausgabe und Speichern
fprintf('\n--- Modell3: Parameteridentifikation abgeschlossen ---\n');
fprintf('Finaler WLS-Fehler: %.4f\n', fval);
for i = 1:numel(p_opt)
    fprintf('%-10s = %.4f\n', namen{i}, p_opt(i));
end

scriptDir   = pwd;
saveDir = fullfile(scriptDir, '..', 'Daten');
save(fullfile(saveDir, 'p_opt_Modell3_woEthO2.mat'), 'p_opt');

%% 4. Visualisierung
t_start = u(1,1);
t_end   = max(Data.Biomasse.t(:)) + 1;
t_sim   = linspace(t_start, t_end, 300);

DOTstern = max(Data.O2.y);
options_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[~, X3] = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p_opt,DOTstern), t_sim, x0, options_ode);

V    = X3(:,1);
cX   = X3(:,2)./V;   cGlc = X3(:,3)./V;   cAm = X3(:,4)./V;
cPh  = X3(:,5)./V;   mB   = X3(:,6);      DOT = X3(:,7);   cEt = X3(:,8)./V;

figure('Name','Modell3 - Parameter Fitting','Position',[200 40 950 1000]);
plot_row(1, Data.Biomasse, t_sim, cX,   'Biomasse',  'c_X (g/L)');
plot_row(2, Data.Glucose,  t_sim, cGlc, 'Glucose',   'c_{Glc} (g/L)');
plot_row(3, Data.Ammonium, t_sim, cAm,  'Ammonium',  'c_{Am} (g/L)');
plot_row(4, Data.Phosphat, t_sim, cPh,  'Phosphat',  'c_{Ph} (g/L)');
plot_row(5, Data.Base,     t_sim, mB,   'Base',      'm_B');
plot_row(6, Data.O2,       t_sim, DOT,  'DOT',       'DOT (%)');
xlabel('BatchAge (h)');



function J = wls_error_m3(p, x0, u, Data)
% WLS über alle 7 Messgrößen. Eine Simulation über das vereinigte
% Zeitraster, danach Zuordnung zu den jeweiligen Messzeitpunkten.

    % Messgrößen sammeln: {t, y, var, Zustandsindex, Divisor}
    % Divisor: true -> Konzentration = mass/V, false -> direkter Zustand
    M = { Data.Biomasse, 2, true;   ...
          Data.Glucose,  3, true;   ...
          Data.Ammonium, 4, true;   ...
          Data.Phosphat, 5, true;   ...
          Data.Base,     6, false;  ...   % mB direkt
          %Data.O2,       7, false;  ...   % DOT direkt
          };

    DOTstern = max(Data.O2.y);

    % Vereinigtes Zeitraster (inkl. Startzeit)
    t0 = u(1,1);
    t_all = t0;
    for i = 1:size(M,1)
        t_all = [t_all; M{i,1}.t(:)];
    end
    t_all = unique(t_all);
    if t_all(1) > t0, t_all = [t0; t_all]; end

    % Einmalige Simulation
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    try
        [~, X] = ode15s(@(t,x) Modell3_woEtOH(t, x, u, p, DOTstern), t_all, x0, opts);
    catch
        J = 1e8; return;
    end
    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:)))
        J = 1e8; return;
    end

    V = X(:,1);
    J = 0;
    for i = 1:size(M,1)
        mess = M{i,1};  idxState = M{i,2};  divByV = M{i,3};
        [~, iT] = ismember(mess.t(:), t_all);

        y_sim = X(iT, idxState);
        if divByV
            y_sim = y_sim ./ V(iT);
        end
        J = J + sum(((mess.y(:) - y_sim).^2) ./ mess.var(:));
    end

    if ~isfinite(J), J = 1e8; end
end


function plot_row(row, mess, t_sim, y_sim, name, ylab)
    subplot(7,1,row);
    errorbar(mess.t, mess.y, sqrt(mess.var), 'o', ...
             'MarkerFaceColor','b','MarkerSize',4); hold on;
    plot(t_sim, y_sim, 'LineWidth', 2);
    title(['Modell3 - ' name]); ylabel(ylab);
    legend('Messung \pm \sigma','Simulation'); grid on;
end