% main_04_parameter_fitting.m
clear; clc; close all;

% 1. Load Preprocessed Data
projectRoot = pwd;
load(fullfile(projectRoot,'..','Daten/Daten_Processed/Processed_FedBatch_Modell3.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
rehash; clear Modell3

Data = TrainData;
u    = Data.u;
x0   = Data.x0;                % [V; mX; mGlc; mAm; mPh; mB; DOT; mEt]


% 1. Anfangswerte und Parameter

%Ausgangslage
namen = {'mumax','KS','YXS', 'YAmX','YPhX','YB_Am','mumax_EtP','mumax_EtX','YGlc_Et','YEt_X','KEt','KGlc_Et','KLa','YXO'};
p0 =    [0.30,   0.50, 0.15,  0.05,  0.02,  1.0,    0.10,        0.10,      0.50,     0.50,   0.10, 0.10,     200,  1.0];
pLB =   [0.01,   0.01, 0.05, 0.001, 0.001,  0.1,   0.001,       0.001,      0.01,     0.01,   0.01, 0.01,     10,   0.1];
pUB =   [1.00,   5.00, 1.00, 1.000, 1.000,  5.0,   1.000,       1.000,      2.00,     2.00,   5.00, 5.00,     800,  5.0];


% YB_Am*YAmX -> YBA
% namen = {'mumax','KS','YXS', 'YPhX','YBA','mumax_EtP','mumax_EtX','YGlc_Et','YEt_X','KEt','KGlc_Et','KLa','YXO'};
% p0 =    [0.30,   0.50, 0.15,  0.02, 0.05,    0.10,        0.10,      0.50,     0.50,   0.10, 0.10,     200,  1.0];
% pLB =   [0.01,   0.01, 0.05, 0.001, 0.0001,   0.001,       0.001,      0.01,     0.01,   0.01, 0.01,     10,   0.1];
% pUB =   [1.00,   5.00, 1.00, 1.000,  5.0,   1.000,       1.000,      2.00,     2.00,   5.00, 5.00,     800,  5.0];





%Optimierung der Startwerte
% Nicht alle startwerte optimieren (V gegeben, DOT0 gemessen, mB = 0)
ix0 = [2 3 4 5];          % mX0, mGlc0, mAm0, mPh0 mitschätzen

                          %  -> V, mB, DOT, mEt bleiben fix (aus Daten)
x0_fit0 = x0(ix0);
x0_fit0 = x0_fit0(:).'; 

%Normierte Parameter
theta0  = ones(size(p0));
thetaLB = pLB ./ p0;
thetaUB = pUB ./ p0;

% --- x0-Anteil normieren (Faktor 0.5 ... 2 um den Startwert) ---
xi0  = ones(1, numel(x0_fit0));   % Start = 1
xiLB = 0.5 * ones(1, numel(x0_fit0));
xiUB = 2.0 * ones(1, numel(x0_fit0));

% --- zusammenfügen: z = [theta, xi] ---
z0 = [theta0,  xi0];
zLB = [thetaLB, xiLB];
zUB = [thetaUB, xiUB];

nP = numel(p0); 

% 2. Parameteridentifikation (WLS)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', 5000, ...
                       'FiniteDifferenceStepSize', 1e-6, ...
                       'OptimalityTolerance', 1e-8, ...
                       'StepTolerance', 1e-10);

obj_fun = @(z) wls_error_m3( z(1:nP).*p0, ...
                             setX0(x0, ix0, z(nP+1:end).*x0_fit0), ...
                             u, Data );

[z_opt, fval] = fmincon(obj_fun, z0, [], [], [], [], ...
                        zLB, zUB, [], options);

% --- Ergebnisse aufteilen ---
p_opt  = z_opt(1:nP) .* p0;
x0_opt = setX0(x0, ix0, z_opt(nP+1:end) .* x0_fit0);


% 3. Ausgabe und Speichern
fprintf('\n--- Modell3: Parameteridentifikation abgeschlossen ---\n');
fprintf('Finaler WLS-Fehler: %.4f\n', fval);
for i = 1:numel(p_opt)
    fprintf('%-10s = %.4f\n', namen{i}, p_opt(i));
end
fprintf('\nOptimierte Anfangswerte:\n');
zn = {'mX0','mGlc0','mAm0','mPh0'};
for k = 1:numel(ix0)
    fprintf('%-6s = %.4f  (Start %.4f)\n', zn{k}, x0_opt(ix0(k)), x0(ix0(k)));
end

scriptDir   = pwd;
saveDir = fullfile(scriptDir, '..', 'Daten');
save(fullfile(saveDir, 'p_opt_Modell3.mat'), 'p_opt');

%% 4. Visualisierung
t_start = u(1,1);
t_end   = max(Data.Biomasse.t(:)) + 1;
t_sim   = linspace(t_start, t_end, 300);

DOTstern = max(Data.O2.y);
options_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[~, X3] = ode15s(@(t,x) Modell3(t,x,u,p_opt,DOTstern), t_sim, x0_opt, options_ode);
%[~, X3] = ode15s(@(t,x) Modell3_YBA(t,x,u,p_opt,DOTstern), t_sim, x0_opt, options_ode);

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
plot_row(7, Data.EtOH,     t_sim, cEt,  'Ethanol',   'c_{Et} (g/L)');
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
          Data.EtOH,     8, true    };

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
        [~, X] = ode15s(@(t,x) Modell3(t, x, u, p, DOTstern), t_all, x0, opts);
        %[~, X] = ode15s(@(t,x) Modell3_YBA(t, x, u, p, DOTstern), t_all, x0, opts);
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
function x0new = setX0(x0base, idx, vals)
    x0new      = x0base;
    x0new(idx) = vals;
end