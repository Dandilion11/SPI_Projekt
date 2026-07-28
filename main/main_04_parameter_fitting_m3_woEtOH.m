% main_04_parameter_fitting.m
clear; clc; close all;

% 1. Load Preprocessed Data
projectRoot = pwd;
load(fullfile(projectRoot,'..','Daten/Daten_Processed/Processed_FedBatch_Modell3.mat'));
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
rehash; clear Modell3_woEtOH

Data = TrainData;
u    = Data.u;
x0   = Data.x0;                % [V; mX; mGlc; mAm; mPh; mB; DOT]


% 1. Anfangswerte und Parameter
%Ausgangslage
namen = {'mumax','KS','YXS', 'YAmX','YPhX','YB_Am','KLa','YXO'};
p0 =    [0.30,   0.50, 0.15,  0.05,  0.02,  1.0,     200,  1.0]; 
pLB =   [0.01,   0.01, 0.05, 0.001, 0.001,  0.1,     10,   0.1];
pUB =   [1.00,   5.00, 1.00, 1.000, 1.000,  5.0,     800,  5.0]; %KS 500? Wert Nimmt Wert von 15 an-> Viel




% 2. Parameteridentifikation (WLS)
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
                       'MaxFunctionEvaluations', 5000, ...
                       'FiniteDifferenceType', 'central', ...       %central umgestellt -> längeres Fitting
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
save(fullfile(saveDir, 'p_opt_Modell3_woEtOH.mat'), 'p_opt');

%% 4. Visualisierung
t_start = u(1,1);
t_end   = max(Data.Biomasse.t(:)) + 1;
t_sim   = linspace(t_start, t_end, 300);

DOTstern = max(Data.O2.y);
options_ode = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[~, X3] = ode15s(@(t,x) Modell3_woEtOH(t,x,u,p_opt,DOTstern), t_sim, x0, options_ode);

V    = X3(:,1);
cX   = X3(:,2)./V;   cGlc = X3(:,3)./V;   cAm = X3(:,4)./V;
cPh  = X3(:,5)./V;   mB   = X3(:,6);      DOT = X3(:,7);

figure('Name','Modell3 - Parameter Fitting','Position',[200 40 950 1000]);
plot_row(1, Data.Biomasse, t_sim, cX,   'Biomasse',  'c_X (g/L)');
plot_row(2, Data.Glucose,  t_sim, cGlc, 'Glucose',   'c_{Glc} (g/L)');
plot_row(3, Data.Ammonium, t_sim, cAm,  'Ammonium',  'c_{Am} (g/L)');
plot_row(4, Data.Phosphat, t_sim, cPh,  'Phosphat',  'c_{Ph} (g/L)');
plot_row(5, Data.Base,     t_sim, mB,   'Base',      'm_B');
plot_row(6, Data.O2,       t_sim, DOT,  'DOT',       'DOT (%)');
xlabel('BatchAge (h)');


function J = wls_error_m3(p, x0, u, Data)
% WLS fuer Modell3 ohne Ethanol.
% Zustaende:
% x = [V; mX; mGlc; mAm; mPh; mB; DOT]
%
% Messgroessen:
% Biomasse  -> mX/V
% Glucose   -> mGlc/V
% Ammonium  -> mAm/V
% Phosphat  -> mPh/V
% Base      -> mB
% DOT       -> DOT

    M = { Data.Biomasse, 'Biomasse',  2, true;   ...
          Data.Glucose,  'Glucose',   3, true;   ...
          Data.Ammonium, 'Ammonium',  4, true;   ...
          Data.Phosphat, 'Phosphat',  5, true;   ...
          Data.Base,     'Base',      6, false;  ...
          Data.O2,       'DOT',       7, false   };

    % Gewichtungsmodus:
    % 'sum'  = klassische WLS: jeder Messpunkt zaehlt einzeln.
    % 'mean' = jede Messgroesse zaehlt etwa gleich stark.
    % Bei vielen DOT-/Online-Punkten ist 'mean' oft stabiler.
    wmode = 'mean';

    % Zusatzgewichtung pro Signal:
    % Reihenfolge: Biomasse, Glucose, Ammonium, Phosphat, Base, DOT
    wsig = [1, 1, 1, 1, 1, 1];

    DOTstern = max(Data.O2.y);

    % Vereinigtes Zeitraster aller Messpunkte
    t_all = [];
    for i = 1:size(M,1)
        t_all = [t_all; M{i,1}.t(:)];
    end

    % Feed-Sprungzeiten zusaetzlich aufnehmen
    t_min = min(t_all);
    t_max = max(t_all);
    tu = u(1,:).';
    tu = tu(tu >= t_min & tu <= t_max);

    t_all = unique([t_all; tu]);

    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

    try
        [~, X] = ode15s(@(t,x) Modell3_woEtOH(t, x, u, p, DOTstern), ...
                        t_all, x0, opts);
    catch
        J = 1e8;
        return;
    end

    if size(X,1) ~= numel(t_all) || any(~isfinite(X(:)))
        J = 1e8;
        return;
    end

    V = X(:,1);
    J = 0;

    for i = 1:size(M,1)
        mess     = M{i,1};
        name     = M{i,2};
        idxState = M{i,3};
        divByV   = M{i,4};

        if isempty(mess.t)
            continue;
        end

        [tf, iT] = ismember(mess.t(:), t_all);
        if any(~tf)
            warning('Nicht alle Messzeitpunkte fuer %s gefunden.', name);
            J = 1e8;
            return;
        end

        y_sim = X(iT, idxState);

        if divByV
            y_sim = y_sim ./ V(iT);
        end

        var_i = max(mess.var(:), eps);

        % Korrekte WLS-Residuen:
        % mess.var ist bereits sigma^2
        r = (mess.y(:) - y_sim) ./ sqrt(var_i);

        switch wmode
            case 'sum'
                contrib = sum(r.^2);

            case 'mean'
                contrib = mean(r.^2);

            otherwise
                error('Unbekannter wmode: %s', wmode);
        end

        J = J + wsig(i) * contrib;
    end

    if ~isfinite(J)
        J = 1e8;
    end
end