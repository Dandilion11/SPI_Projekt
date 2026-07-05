%Loads .mat files, filters batch vs. fed-batch phases, and formats u matrices
clear; clc; close all;

% 1. Pfade definieren (dynamisch oder relativ anpassen)
scriptDir = pwd;
pfad_03 = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef03.mat');
pfad_04 = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef04.mat');


% 2. Daten laden
daten_03 = load(pfad_03);
daten_04 = load(pfad_04);

exp_03 = daten_03.Mess;
exp_04 = daten_04.Mess;


% 1. Koeffizienten für Varianzberechnung definieren -> GEGENCHECKEN!!!
a_bio = 0.02; b_bio = 0.015;
a_glc = 0.06; b_glc = 0.25;

a_o2 = 0.02; 
b_o2 = 0.5;

a_am   = 0.05; b_am   = 0.02;
a_base = 0.05; b_base = 0.01;

% 2. Funktion zum Filtern und Varianz berechnen (t <= 20h)
filterAndVarBatch = @(t, y, a, b) struct('t', t(t <= 20), ...
    'y', y(t <= 20), ...
    'var', (a .* y(t <= 20) + b).^2);

% 3. Trainingsdaten extrahieren und Varianz berechnen (RamScDef03)
TrainData.Biomasse = filterAndVarBatch(exp_03.Messdaten.Biomasse.BatchAge, ...
    exp_03.Messdaten.Biomasse.Wert, ...
    a_bio, b_bio);

TrainData.Glucose = filterAndVarBatch(exp_03.Messdaten.Glucose.BatchAge, ...
    exp_03.Messdaten.Glucose.Wert, ...
    a_glc, b_glc);

TrainData.O2 = filterAndVarBatch(exp_03.Messdaten.O2.BatchAge, ...
    exp_03.Messdaten.O2.Wert, ...
    a_o2, b_o2);

TrainData.Ammonium = filterAndVarBatch(exp_03.Messdaten.Ammonium.BatchAge, ...
    exp_03.Messdaten.Ammonium.Wert, a_am, b_am);

TrainData.Base = filterAndVarBatch(exp_03.Messdaten.BASE.BatchAge, ...
    exp_03.Messdaten.BASE.Wert, a_base, b_base);

% 4. Validierungsdaten extrahieren und Varianz berechnen (RamScDef04)
ValData.Biomasse = filterAndVarBatch(exp_04.Messdaten.Biomasse.BatchAge, ...
    exp_04.Messdaten.Biomasse.Wert, ...
    a_bio, b_bio);

ValData.Glucose = filterAndVarBatch(exp_04.Messdaten.Glucose.BatchAge, ...
    exp_04.Messdaten.Glucose.Wert, ...
    a_glc, b_glc);

ValData.O2 = filterAndVarBatch(exp_04.Messdaten.O2.BatchAge, ...
    exp_04.Messdaten.O2.Wert, ...
    a_o2, b_o2);

ValData.Ammonium = filterAndVarBatch(exp_04.Messdaten.Ammonium.BatchAge, ...
    exp_04.Messdaten.Ammonium.Wert, a_am, b_am);

ValData.Base = filterAndVarBatch(exp_04.Messdaten.BASE.BatchAge, ...
    exp_04.Messdaten.BASE.Wert, a_base, b_base);

% Dichte p = 1 kg/L wird angenommen. Daher: Volumen (L) = Gewicht (kg)
TrainData.V = exp_03.Messdaten.Weight.Wert(1);
ValData.V   = exp_04.Messdaten.Weight.Wert(1);

% 1. Pfad für den neuen Ordner definieren (eine Ebene über dem Skript-Ordner)
saveDir = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed');

% 2. Überprüfen, ob der Ordner existiert, falls nicht: erstellen
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% 3. Vollständigen Dateipfad für die Speicherung definieren
savePath = fullfile(saveDir, 'Processed_Batch_Data.mat');

% 4. Daten speichern
save(savePath, 'TrainData', 'ValData');
fprintf('Daten erfolgreich gespeichert in:\n%s\n', savePath);



disp('Trainings- und Validierungsdaten erfolgreich geladen und auf t <= 20h beschnitten.');