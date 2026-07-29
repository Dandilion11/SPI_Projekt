%% main_01_data_preprocessing_m3.m
%  Preprocessing Modell 3 (Fed-Batch, ohne Ethanol) mit MEHREREN
%  Trainingsdatensaetzen und EINEM Validierungsdatensatz.
%
%  Nutzbare Fed-Batch-Experimente: 07, 10, 06.
%  Default: Training = {10, 07}, Validierung = 06.
%
%  Ausgabe (Processed_FedBatch_Modell3.mat):
%    TrainData : 1xN Cell-Array, ein Struct je Trainingsexperiment
%    ValData   : ein Struct (Validierungsexperiment)
%  Jeder Struct: .Name .Biomasse .Glucose .Ammonium .Phosphat .O2 .Base
%                .V .u .x0 .DOTstern
% =========================================================================
clear; clc; close all;

%% === Konfiguration ===
trainNums = [10 7];     % Trainingsexperimente
valNum    = 6;          % Validierungsexperiment

filepath    = pwd;
datenordner = fullfile(filepath, '..', 'Daten', 'MessDaten_SPI1_Projekt');

% a/b aus dem Kraemer und King 2017 Paper
ab.Biomasse = [0.02  0.015];
ab.Glucose  = [0.06  0.25 ];
ab.Ammonium = [0.06  0.01 ];
ab.Phosphat = [0.07  0.01 ];
ab.O2       = [0.02  0.5  ]; % Muss nochmal recherchiert werden (nicht in Kraemer/King)
ab.Base     = [0.01  10   ];

%% Trainingsdaten (Cell-Array: ein Eintrag pro Experiment)
TrainData = cell(1, numel(trainNums));
for i = 1:numel(trainNums)
    TrainData{i} = aufbereiten(trainNums(i), datenordner, ab);
    fprintf('Training  %s geladen (Biomasse: %d Punkte).\n', ...
            TrainData{i}.Name, numel(TrainData{i}.Biomasse.t));
end

%% Validierungsdaten (ein Experiment)
ValData = aufbereiten(valNum, datenordner, ab);
fprintf('Validierung %s geladen (Biomasse: %d Punkte).\n', ...
        ValData.Name, numel(ValData.Biomasse.t));

%% Speichern
zielordner = fullfile(filepath, '..', 'Daten', 'Daten_Processed');
if ~exist(zielordner, 'dir'); mkdir(zielordner); end
ziel = fullfile(zielordner, 'Processed_FedBatch_Modell3.mat');
save(ziel, 'TrainData', 'ValData', 'trainNums', 'valNum');

fprintf('\nFertig. Gespeichert unter:\n%s\n', ziel);
fprintf('Training = %s | Validierung = %d\n', mat2str(trainNums), valNum);


%% ========================================================================
%  Lokale Funktionen
%% ========================================================================
function D = aufbereiten(num, ordner, ab)
    fname = fullfile(ordner, sprintf('Mess_RamScDef%02d.mat', num));
    Mess  = load(fname).Mess;
    M     = Mess.Messdaten;

    D.Name     = sprintf('RamScDef%02d', num);
    D.Biomasse = messwert(M.Biomasse, ab.Biomasse);
    D.Glucose  = messwert(M.Glucose,  ab.Glucose);
    D.Ammonium = messwert(M.Ammonium, ab.Ammonium);
    D.Phosphat = messwert(M.Phosphat, ab.Phosphat);
    D.O2       = messwert(M.O2,        ab.O2);     % = DOT in %
    D.Base     = messwert(M.BASE,      ab.Base);

    % Volumen (Dichte 1 kg/L angenommen -> V[L] = Gewicht[kg])
    D.V.t = M.Weight.BatchAge(:);
    D.V.y = M.Weight.Wert(:);
    V0    = D.V.y(1);

    % Stellgroessen-Matrix (Zeile1 = BatchAge, Feeds, Base, Saeure)
    D.u = Mess.u;

    % DOTstern (aus vollem O2-Signal), robust gegen spaeteres Fenstern
    D.DOTstern = max(D.O2.y);

    % Startwerte fuer die Simulation: [V; mX; mGlc; mAm; mPh; mB; DOT]
    D.x0 = [ V0;                       % V
             D.Biomasse.y(1) * V0;     % mX
             D.Glucose.y(1)  * V0;     % mGlc
             D.Ammonium.y(1) * V0;     % mAm
             D.Phosphat.y(1) * V0;     % mPh
             0;                        % mB (kumulierte Base startet bei 0)
             D.O2.y(1) ];              % DOT
end


function s = messwert(feld, ab)
% Zeit + Wert raus, Messvarianz (a*y+b)^2
    t = feld.BatchAge(:);
    y = feld.Wert(:);
    s.t   = t;
    s.y   = y;
    s.var = (ab(1).*y + ab(2)).^2;
end
