%Loads .mat files, filters batch vs. fed-batch phases, and formats u matrices
clear; clc; close all;

% === Schalter: Varianzquelle ===
% true  -> Varianz aus dem Messdaten-Struct (feld.Variance)
% false -> Varianz berechnet ueber (a*y + b)^2
useMeasuredVar = true;

% 1. Pfade definieren (dynamisch oder relativ anpassen)
scriptDir = pwd;
pfad_03 = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef03.mat');
pfad_04 = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt', 'Mess_RamScDef04.mat');

% 2. Daten laden
daten_03 = load(pfad_03);
daten_04 = load(pfad_04);

exp_03 = daten_03.Mess;
exp_04 = daten_04.Mess;

% 3. Koeffizienten fuer Varianzberechnung definieren -> GEGENCHECKEN!!!
%    (a, b)
ab.Biomasse = [0.02  0.015];
ab.Glucose  = [0.06  0.25 ];
ab.O2       = [0.02  0.5  ];
ab.Ammonium = [0.05  0.02 ];
ab.Base     = [0.05  0.01 ];

% 4. Trainingsdaten extrahieren und Varianz berechnen (RamScDef03)
TrainData.Biomasse = messwert(exp_03.Messdaten.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar);
TrainData.Glucose  = messwert(exp_03.Messdaten.Glucose,  ab.Glucose,  "Glucose",  useMeasuredVar);
TrainData.O2       = messwert(exp_03.Messdaten.O2,        ab.O2,       "O2",       useMeasuredVar);
TrainData.Ammonium = messwert(exp_03.Messdaten.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar);
TrainData.Base     = messwert(exp_03.Messdaten.BASE,      ab.Base,     "Base",     useMeasuredVar);

% 5. Validierungsdaten extrahieren und Varianz berechnen (RamScDef04)
ValData.Biomasse = messwert(exp_04.Messdaten.Biomasse, ab.Biomasse, "Biomasse", useMeasuredVar);
ValData.Glucose  = messwert(exp_04.Messdaten.Glucose,  ab.Glucose,  "Glucose",  useMeasuredVar);
ValData.O2       = messwert(exp_04.Messdaten.O2,        ab.O2,       "O2",       useMeasuredVar);
ValData.Ammonium = messwert(exp_04.Messdaten.Ammonium, ab.Ammonium, "Ammonium", useMeasuredVar);
ValData.Base     = messwert(exp_04.Messdaten.BASE,      ab.Base,     "Base",     useMeasuredVar);

TrainData.V = exp_03.Messdaten.Weight.Wert(1);
ValData.V   = exp_04.Messdaten.Weight.Wert(1);

% 6. Speicherordner definieren (eine Ebene ueber dem Skript-Ordner)
saveDir = fullfile(scriptDir, '..', 'Daten', 'Daten_Processed');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% 7. Vollstaendigen Dateipfad definieren und speichern
savePath = fullfile(saveDir, 'Processed_Batch_Data.mat');
save(savePath, 'TrainData', 'ValData');
fprintf('Daten erfolgreich gespeichert in:\n%s\n', savePath);

disp('Trainings- und Validierungsdaten erfolgreich geladen und auf t <= 20h beschnitten.');


%% ---------------------------------------------------------------
function s = messwert(feld, ab, name, useMeasuredVar)
% Filtert auf Batch-Phase (t <= 20h), zieht Zeit + Wert raus
% und liefert die Messvarianz entweder aus dem Struct oder berechnet.

t = feld.BatchAge(:);
y = feld.Wert(:);

% Batch-Phase filtern (t <= 20h)
idx = t <= 20;
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