% main_01_data_preprocessing.m
%
% Preprocessing fuer die Batch-Modelle 1 und 2.
%   Training:    RamScDef03, Batch-Fenster t <= 6 h
%   Validierung: RamScDef04, Batch-Fenster t <= 8 h
%
% Nach diesen Zeiten wird in beiden Experimenten Glucose zugefuettert.
% Modell 1 und 2 haben keine Eingaenge, deshalb wird nur die feedfreie
% Batch-Phase verwendet.

clear; clc; close all;

% Varianzquelle: false -> aus den a/b-Parametern des Datensatzes berechnen,
% true -> die mitgelieferte Variance-Zeitreihe verwenden.
useMeasuredVar = false;

tCut_03 = 6;    % Batch-Fenster Training   [h]
tCut_04 = 8;    % Batch-Fenster Validierung [h]

%% Laden ------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));

datenordner = fullfile(scriptDir, '..', '..', 'Daten', 'MessDaten_SPI1_Projekt');
exp_03 = load(fullfile(datenordner,'Mess_RamScDef03.mat')).Mess;
exp_04 = load(fullfile(datenordner,'Mess_RamScDef04.mat')).Mess;

%% Aufbereiten ------------------------------------------------------------
% Kanalname im Ergebnis -> Feldname im Mess-Struct
kanaele = {'Biomasse','Biomasse'; 'Glucose','Glucose'; 'O2','O2'; ...
           'Ammonium','Ammonium'; 'Base','BASE'};

for i = 1:size(kanaele,1)
    aus = kanaele{i,1};  ein = kanaele{i,2};
    TrainData.(aus) = messwert(exp_03.Messdaten.(ein), aus, useMeasuredVar, tCut_03);
    ValData.(aus)   = messwert(exp_04.Messdaten.(ein), aus, useMeasuredVar, tCut_04);
end

% Startvolumen (Dichte 1 kg/L angenommen -> V[L] = Gewicht[kg])
TrainData.V = exp_03.Messdaten.Weight.Wert(1);
ValData.V   = exp_04.Messdaten.Weight.Wert(1);

%% Speichern --------------------------------------------------------------
saveDir = fullfile(scriptDir, '..', '..', 'Daten','Daten_Processed');
if ~exist(saveDir,'dir'); mkdir(saveDir); end
savePath = fullfile(saveDir,'Processed_Batch_Data.mat');
save(savePath, 'TrainData', 'ValData');

fprintf('Gespeichert: %s\n', savePath);
fprintf('Training   RamScDef03 (t <= %.0f h): %d Biomasse-Punkte\n', ...
        tCut_03, numel(TrainData.Biomasse.t));
fprintf('Validierung RamScDef04 (t <= %.0f h): %d Biomasse-Punkte\n', ...
        tCut_04, numel(ValData.Biomasse.t));
fprintf('sigma: Biomasse %.4f | Ammonium %.4f | Base %.5f\n', ...
        median(sqrt(TrainData.Biomasse.var)), ...
        median(sqrt(TrainData.Ammonium.var)), ...
        median(sqrt(TrainData.Base.var)));


%% ======================================================================
%  Hilfsfunktion
%% ======================================================================

function s = messwert(feld, name, useMeasuredVar, tMax)
% Schneidet auf die Batch-Phase, raeumt auf und berechnet die Messvarianz.
%
% Sortieren ist zwingend: die Fitting-Skripte nehmen x0 = y(1), also den
% ersten Wert NACH der Zeit -- nicht den ersten im Speicher.
%
% Varianzmodell: sigma^2 = a*y + b mit a,b aus dem Datensatz. a und b sind
% in den Doku-Einheiten kalibriert, die Werte liegen hier aber schon in
% Modell-Einheiten vor. Varianz skaliert QUADRATISCH mit einer
% Einheitenumrechnung, daher /1e6 (mL->L) bzw. *1e4 (Anteil->%).

    t = feld.BatchAge(:);
    y = feld.Wert(:);

    v = nan(size(y));
    if isfield(feld,'Variance') && ~isempty(feld.Variance)
        n = min(numel(v), numel(feld.Variance));
        v(1:n) = feld.Variance(1:n);
    end

    % Die Felder koennen unterschiedlich lang sein -> gemeinsamer Bereich
    n = min([numel(t) numel(y) numel(v)]);
    t = t(1:n);  y = y(1:n);  v = v(1:n);

    ok = isfinite(t) & isfinite(y) & t <= tMax;
    t = t(ok);  y = y(ok);  v = v(ok);
    [t, ord] = sort(t);  y = y(ord);  v = v(ord);

    s.t = t;
    s.y = y;

    if useMeasuredVar
        switch name
            case "Base", s.var = v / 1e6;      % mL -> L
            %case "O2",   s.var = v * 1e4;      % Anteil -> %
            otherwise,   s.var = v;
        end
    else
        a = feld.VarParam.a;  b = feld.VarParam.b;
        switch name
            case "Base", s.var = a.*y/1000 + b/1e6;
            %case "O2",   s.var = 100*a.*y + 1e4*b;
            otherwise,   s.var = a.*y + b;
        end
    end
    s.var = max(s.var, eps);
end