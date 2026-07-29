%% main_01_data_preprocessing.m
%  Preprocessing fuer die Batch-Modelle 1 & 2 mit MEHREREN Trainings-
%  datensaetzen und EINEM Validierungsdatensatz.
%
%  Nutzbare Batch-Experimente: 03, 04, 06   (05 raus -> keine O2-Messung).
%  Default: Training = {03, 04}, Validierung = 06.
%
%  Ausgabe (Processed_Batch_Data.mat):
%    TrainData : 1xN Cell-Array, ein Struct je Trainingsexperiment
%    ValData   : ein Struct (Validierungsexperiment)
%  Jeder Struct: .Name .Biomasse .Glucose .O2 .Ammonium .Base .V
%  (jede Messgroesse mit Feldern .t .y .var)
%
%  - Batch-Phase bei t <= tMax (10 h) abgeschnitten (danach nur Online-Punkte).
%  - Exp03: Glucose steigt ab ~6 h ohne Fuetterung (Messartefakt) -> ab
%    tGlcSplit exponentiell extrapoliert.
% =========================================================================
clear; clc; close all;

%% === Konfiguration ===
trainNums = [3 4];      % Trainingsexperimente (Nummern der RamScDef-Dateien)
valNum    = 6;          % Validierungsexperiment

useMeasuredVar = false; % true: feld.Variance | false: (a*y+b)^2
tMax          = 10;     % [h] Batch-Grenze
glcExtrapNums = 3;      % Experimente, deren Glucose extrapoliert wird
tGlcSplit     = 6;      % [h] ab hier Glucose exponentiell extrapolieren

%% Pfade
scriptDir   = pwd;
datenordner = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt');

%% Koeffizienten fuer Varianzberechnung (a, b) -> GEGENCHECKEN!!!
ab.Biomasse = [0.02  0.015];
ab.Glucose  = [0.06  0.25 ];
ab.O2       = [0.02  0.5  ];
ab.Ammonium = [0.05  0.02 ];
ab.Base     = [0.05  0.01 ];

%% Trainingsdaten (Cell-Array: ein Eintrag pro Experiment)
TrainData = cell(1, numel(trainNums));
for i = 1:numel(trainNums)
    n = trainNums(i);
    TrainData{i} = aufbereiten(n, datenordner, ab, useMeasuredVar, tMax, ...
                               ismember(n, glcExtrapNums), tGlcSplit);
    fprintf('Training  %s geladen (Biomasse: %d Punkte).\n', ...
            TrainData{i}.Name, numel(TrainData{i}.Biomasse.t));
end

%% Validierungsdaten (ein Experiment)
ValData = aufbereiten(valNum, datenordner, ab, useMeasuredVar, tMax, ...
                      ismember(valNum, glcExtrapNums), tGlcSplit);
fprintf('Validierung %s geladen (Biomasse: %d Punkte).\n', ...
        ValData.Name, numel(ValData.Biomasse.t));

%% Speichern
saveDir = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end
savePath = fullfile(saveDir, 'Processed_Batch_Data.mat');
save(savePath, 'TrainData', 'ValData', 'trainNums', 'valNum');

fprintf('\nGespeichert in:\n%s\n', savePath);
fprintf('Training = %s | Validierung = %d | Cut t <= %.0f h\n', ...
        mat2str(trainNums), valNum, tMax);


%% ========================================================================
%  Lokale Funktionen
%% ========================================================================
function D = aufbereiten(num, ordner, ab, useMeasuredVar, tMax, doGlcExtrap, tGlcSplit)
% Laedt ein Experiment und extrahiert alle Messgroessen fuer Modell 1 & 2.
    fname = fullfile(ordner, sprintf('Mess_RamScDef%02d.mat', num));
    Mess  = load(fname).Mess;
    M     = Mess.Messdaten;

    D.Name     = sprintf('RamScDef%02d', num);
    D.Biomasse = messwert(M.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar, tMax);
    if doGlcExtrap
        D.Glucose = messwert_glc(M.Glucose, ab.Glucose, useMeasuredVar, tMax, tGlcSplit);
    else
        D.Glucose = messwert(M.Glucose, ab.Glucose, "Glucose", useMeasuredVar, tMax);
    end
    D.O2       = messwert(M.O2,       ab.O2,       "O2",       useMeasuredVar, tMax);
    D.Ammonium = messwert(M.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar, tMax);
    D.Base     = messwert(M.BASE,     ab.Base,     "Base",     useMeasuredVar, tMax);
    D.V        = M.Weight.Wert(1);
end


function s = messwert(feld, ab, name, useMeasuredVar, tMax)
% Filtert auf Batch-Phase (t <= tMax) und liefert die Messvarianz.
    t = feld.BatchAge(:);
    y = feld.Wert(:);
    idx = t <= tMax;
    s.t = t(idx);
    s.y = y(idx);
    if useMeasuredVar
        v = feld.Variance(:);  v = v(idx);
        if     name == "Base", v = v / 1000;   % mL -> L
        elseif name == "O2",   v = v * 100;    % -> %
        end
        s.var = v;
    else
        s.var = (ab(1) .* s.y + ab(2)).^2;
    end
end


function s = messwert_glc(feld, ab, useMeasuredVar, tMax, tSplit)
% Wie messwert(), aber Glucose ab t > tSplit exponentiell extrapoliert:
% Modell y = A*exp(k*t), log-linear gefittet ueber die Batch-Phase t <= tSplit.
    t = feld.BatchAge(:);
    y = feld.Wert(:);
    idx = t <= tMax;  t = t(idx);  y = y(idx);
    if useMeasuredVar
        v = feld.Variance(:);  v = v(idx);
    else
        v = (ab(1) .* y + ab(2)).^2;
    end

    mask = t <= tSplit;
    tb = t(mask);  yb = y(mask);  pos = yb > 0;
    if nnz(pos) >= 2
        c   = polyfit(tb(pos), log(yb(pos)), 1);   % c(1)=k, c(2)=log(A)
        ext = ~mask;
        y(ext) = max(exp(polyval(c, t(ext))), 0);  % y = A*exp(k*t) >= 0
        v(ext) = (ab(1) .* y(ext) + ab(2)).^2;
    else
        warning('Glucose: zu wenige positive Punkte (t <= %.1f h) fuer exp. Extrapolation.', tSplit);
    end

    s.t   = t;
    s.y   = y;
    s.var = v;
end
