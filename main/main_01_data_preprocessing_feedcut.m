%% main_01_data_preprocessing_feedcut.m
%  Preprocessing fuer die Batch-Modelle 1 & 2 (analog zu
%  main_01_data_preprocessing.m).
%
%  - Batch-Phase bei t <= 10 h abgeschnitten.
%  - SONDERFALL Exp03 (Training): Die gemessene Glucose STEIGT ab ca. 6 h
%    wieder an, obwohl NICHT gefuettert wird. Diese Werte
%    werden ab tGlcSplit = 6 h verworfen und durch eine Extrapolation aus
%    der Batch-Phase (t <= 6 h) ersetzt (linear, auf >= 0 begrenzt).
%
%  Training:    RamScDef03      Validierung: RamScDef04
% =========================================================================
clear; clc; close all;

% === Schalter: Varianzquelle ===
% true  -> Varianz aus dem Messdaten-Struct (feld.Variance)
% false -> Varianz berechnet ueber (a*y + b)^2
useMeasuredVar = true;

% === Glucose-Extrapolation (nur Exp03) ===
tMax       = 10;     % [h] Batch-Grenze wie im Original
tGlcSplit  = 6;      % [h] ab hier Glucose von Exp03 extrapolieren statt messen

% 1. Pfade definieren
scriptDir = pwd;
pfad_03 = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef03.mat');
pfad_04 = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef04.mat');

% 2. Daten laden
daten_03 = load(pfad_03);
daten_04 = load(pfad_04);

exp_03 = daten_03.Mess;
exp_04 = daten_04.Mess;

% 3. Koeffizienten fuer Varianzberechnung (a, b) -> GEGENCHECKEN!!!
ab.Biomasse = [0.02  0.015];
ab.Glucose  = [0.06  0.25 ];
ab.O2       = [0.02  0.5  ];
ab.Ammonium = [0.05  0.02 ];
ab.Base     = [0.05  0.01 ];

% 4. Trainingsdaten extrahieren und Varianz berechnen (RamScDef03)
TrainData.Biomasse = messwert(exp_03.Messdaten.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar, tMax);
TrainData.Glucose  = messwert_glc(exp_03.Messdaten.Glucose, ab.Glucose, useMeasuredVar, tMax, tGlcSplit);
TrainData.O2       = messwert(exp_03.Messdaten.O2,        ab.O2,       "O2",       useMeasuredVar, tMax);
TrainData.Ammonium = messwert(exp_03.Messdaten.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar, tMax);
TrainData.Base     = messwert(exp_03.Messdaten.BASE,      ab.Base,     "Base",     useMeasuredVar, tMax);

% 5. Validierungsdaten extrahieren und Varianz berechnen (RamScDef04)
ValData.Biomasse = messwert(exp_04.Messdaten.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar, tMax);
ValData.Glucose  = messwert(exp_04.Messdaten.Glucose,  ab.Glucose,  "Glucose",  useMeasuredVar, tMax);
ValData.O2       = messwert(exp_04.Messdaten.O2,        ab.O2,       "O2",       useMeasuredVar, tMax);
ValData.Ammonium = messwert(exp_04.Messdaten.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar, tMax);
ValData.Base     = messwert(exp_04.Messdaten.BASE,      ab.Base,     "Base",     useMeasuredVar, tMax);

TrainData.V = exp_03.Messdaten.Weight.Wert(1);
ValData.V   = exp_04.Messdaten.Weight.Wert(1);

% 6. Speicherordner definieren
saveDir = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% 7. Speichern
savePath = fullfile(saveDir, 'Processed_Batch_Data.mat');
save(savePath, 'TrainData', 'ValData');
fprintf('Daten erfolgreich gespeichert in:\n%s\n', savePath);
disp('Trainings- und Validierungsdaten geladen und auf t <= 10h beschnitten.');
fprintf('Exp03-Glucose ab %.1f h extrapoliert.\n', tGlcSplit);


%% ---------------------------------------------------------------
function s = messwert(feld, ab, name, useMeasuredVar, tMax)
% Filtert auf Batch-Phase (t <= tMax), zieht Zeit + Wert raus
% und liefert die Messvarianz entweder aus dem Struct oder berechnet.

t = feld.BatchAge(:);
y = feld.Wert(:);

% Batch-Phase filtern
idx = t <= tMax;
s.t = t(idx);
s.y = y(idx);

if useMeasuredVar
    % Varianz aus dem Messdaten-Struct
    v = feld.Variance(:);
    v = v(idx);
    if name == "Base"
        v = v / 1000;   % Umrechnung von mL in L
    elseif name == "O2"
        v = v * 100; % Umrechnung in %
    end
    s.var = v;
else
    s.var = (ab(1) .* s.y + ab(2)).^2;
end
end


%% ---------------------------------------------------------------
function s = messwert_glc(feld, ab, useMeasuredVar, tMax, tSplit)
% Wie messwert(), aber fuer die Glucose von Exp03:
% Ab t > tSplit werden die (fehlerhaft ansteigenden) Messwerte durch eine
% lineare Extrapolation der Batch-Phase (t <= tSplit) ersetzt (>= 0).

t = feld.BatchAge(:);
y = feld.Wert(:);

% Batch-Phase filtern (t <= tMax)
idx = t <= tMax;
t = t(idx);
y = y(idx);

if useMeasuredVar
    v = feld.Variance(:);
    v = v(idx);
else
    v = (ab(1) .* y + ab(2)).^2;
end

% Lineare Trendgerade aus der verwertbaren Batch-Phase (t <= tSplit)
mask = t <= tSplit;
if nnz(mask) >= 2
    pfit = polyfit(t(mask), y(mask), 1);   % [Steigung, Achsenabschnitt]
    ext  = ~mask;                          % zu extrapolierende Punkte (t > tSplit)
    y(ext) = max(polyval(pfit, t(ext)), 0);% ersetzen, auf >= 0 begrenzen
    % Varianz der extrapolierten Punkte als Schaetzung ueber (a*y+b)^2
    v(ext) = (ab(1) .* y(ext) + ab(2)).^2;
else
    warning('Exp03-Glucose: zu wenige Punkte (t <= %.1f h) fuer Extrapolation.', tSplit);
end

s.t   = t;
s.y   = y;
s.var = v;
end
