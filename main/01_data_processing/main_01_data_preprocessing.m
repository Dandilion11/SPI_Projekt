%% main_01_data_preprocessing.m
%  Preprocessing fuer die Batch-Modelle 1 & 2.
%  Training:    RamScDef03  -> Batch-Fenster t <= 6 h
%  Validierung: RamScDef04  -> Batch-Fenster t <= 8 h
%
%  Begruendung der Cuts: In beiden Experimenten wird nach diesen Zeiten
%  Glucose zugefuettert (Glucose steigt wieder an). Da Modell 1/2 keine
%  Eingaenge haben, wird nur die feedfreie Batch-Phase verwendet.
% =========================================================================
clear; clc; close all;

% === Schalter: Varianzquelle ===
% true  -> Varianz aus dem Messdaten-Struct (feld.Variance)
% false -> Varianz berechnet ueber (a*y + b)^2
useMeasuredVar = false;

% === Batch-Fenster je Experiment [h] ===
tCut_03 = 6;    % Training  RamScDef03
tCut_04 = 8;    % Validierung RamScDef04

% 1. Pfade
scriptDir = pwd;
pfad_03 = fullfile(scriptDir, 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef03.mat');
pfad_04 = fullfile(scriptDir, 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef04.mat');

% 2. Laden
exp_03 = load(pfad_03).Mess;
exp_04 = load(pfad_04).Mess;

% 3. Koeffizienten fuer Varianzberechnung (a, b) -> GEGENCHECKEN!!!
ab.Biomasse = [0.02  0.015];
ab.Glucose  = [0.06  0.25 ];
ab.O2       = [0.02  0.5  ];
ab.Ammonium = [0.05  0.02 ];
ab.Base     = [0.05  0.01 ];

% 4. Trainingsdaten (RamScDef03, t <= tCut_03)
TrainData.Biomasse = messwert(exp_03.Messdaten.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar, tCut_03);
TrainData.Glucose  = messwert(exp_03.Messdaten.Glucose,  ab.Glucose,  "Glucose",  useMeasuredVar, tCut_03);
TrainData.O2       = messwert(exp_03.Messdaten.O2,        ab.O2,       "O2",       useMeasuredVar, tCut_03);
TrainData.Ammonium = messwert(exp_03.Messdaten.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar, tCut_03);
TrainData.Base     = messwert(exp_03.Messdaten.BASE,      ab.Base,     "Base",     useMeasuredVar, tCut_03);

% 5. Validierungsdaten (RamScDef04, t <= tCut_04)
ValData.Biomasse = messwert(exp_04.Messdaten.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar, tCut_04);
ValData.Glucose  = messwert(exp_04.Messdaten.Glucose,  ab.Glucose,  "Glucose",  useMeasuredVar, tCut_04);
ValData.O2       = messwert(exp_04.Messdaten.O2,        ab.O2,       "O2",       useMeasuredVar, tCut_04);
ValData.Ammonium = messwert(exp_04.Messdaten.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar, tCut_04);
ValData.Base     = messwert(exp_04.Messdaten.BASE,      ab.Base,     "Base",     useMeasuredVar, tCut_04);


TrainData.V = exp_03.Messdaten.Weight.Wert(1);
ValData.V   = exp_04.Messdaten.Weight.Wert(1);

% 6. Speichern
saveDir = fullfile(scriptDir, 'Daten', 'Daten_Processed');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end
savePath = fullfile(saveDir, 'Processed_Batch_Data.mat');
save(savePath, 'TrainData', 'ValData');
fprintf('Daten gespeichert in:\n%s\n', savePath);
fprintf('Training RamScDef03 (t <= %.0f h): %d Biomasse-Punkte\n', tCut_03, numel(TrainData.Biomasse.t));
fprintf('Validierung RamScDef04 (t <= %.0f h): %d Biomasse-Punkte\n', tCut_04, numel(ValData.Biomasse.t));


function s = messwert(feld, ab, name, useMeasuredVar, tMax)
% Filtert auf die Batch-Phase (t <= tMax) und liefert die Messvarianz.
t = feld.BatchAge(:);
y = feld.Wert(:);

idx = t <= tMax;
s.t = t(idx);
s.y = y(idx);

if useMeasuredVar
    v = feld.Variance(:);
    v = v(idx);
    if name == "Base"
        v = v / 1000;   % Umrechnung von mL in L
    elseif name == "O2"
        v = v * 100;    % Umrechnung in %
    end
    s.var = v;
else
    s.var = (ab(1) .* s.y + ab(2)).^2;
end
end

function probe = probenahme(Data, tMax)

    t = Data.Probenahmen.BatchAge;

    idx = t <= tMax;

    probe.t = t(idx);
    probe.V = Data.Probenahmen.Volumen(idx);

end

