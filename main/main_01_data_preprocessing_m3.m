%% Preprocessing Modell 3 (Fed-Batch mit Ethanol)
% Training:    RamScDef10
% Validierung: RamScDef07
clear; clc; close all;

filepath = pwd;
datenordner = fullfile(filepath, '..', 'Daten', 'MessDaten_SPI1_Projekt');

train = load(fullfile(datenordner, 'Mess_RamScDef10.mat')).Mess;   % zum Trainieren
val   = load(fullfile(datenordner, 'Mess_RamScDef07.mat')).Mess;   % zum Validieren

% a/b aus dem Krämer und King 2017 Paper
ab.Biomasse = [0.02  0.015];
ab.Glucose  = [0.06  0.25 ];
ab.Ammonium = [0.06  0.01 ];
ab.Phosphat = [0.07  0.01 ];
ab.O2       = [0.02  0.5  ]; % Muss nochmal recherchiert werden. Steht nicht in Krämer und King
ab.Base     = [0.01  10   ];

% Training und Validierung mit derselben Routine aufbereiten
TrainData = aufbereiten(train, ab);
ValData   = aufbereiten(val,   ab);

% speichern
zielordner = fullfile(filepath, '..', 'Daten', 'Daten_Processed');
if ~exist(zielordner, 'dir'); mkdir(zielordner); end
ziel = fullfile(zielordner, 'Processed_FedBatch_Modell3.mat');
save(ziel, 'TrainData', 'ValData');

fprintf('Fertig. Gespeichert unter:\n%s\n', ziel);


%% ---------------------------------------------------------------
function D = aufbereiten(Mess, ab)

M = Mess.Messdaten;
% Messgroessen einlesen (Zeit in h, Wert in g/L bzw. % bei O2)
D.Biomasse = messwert(M.Biomasse, ab.Biomasse, "Biomasse");
D.Glucose  = messwert(M.Glucose,  ab.Glucose, "Glucose");
D.Ammonium = messwert(M.Ammonium, ab.Ammonium, "Ammonium");
D.Phosphat = messwert(M.Phosphat, ab.Phosphat, "Phosphat");
D.O2       = messwert(M.O2,       ab.O2, "O2");     % = DOT in %
D.Base     = messwert(M.BASE,     ab.Base, "Base");

% Volumen (Dichte 1 kg/L angenommen -> V[L] = Gewicht[kg])
D.V.t = M.Weight.BatchAge(:);
D.V.y = M.Weight.Wert(:);
V0    = D.V.y(1);

% Stellgroessen-Matrix (Zeile1 = BatchAge, Feeds, Base, Saeure)
D.u = Mess.u;

% Startwerte fuer die Simulation zusammenbauen
% Reihenfolge wie in Modell3: [V; mX; mGlc; mAm; mPh; mB; DOT; mEt]
cX0   = D.Biomasse.y(1);
cGlc0 = D.Glucose.y(1);
cAm0  = D.Ammonium.y(1);
cPh0  = D.Phosphat.y(1);
DOT0  = D.O2.y(1);

D.x0 = [ V0;            % V
         cX0   * V0;    % mX
         cGlc0 * V0;    % mGlc
         cAm0  * V0;    % mAm
         cPh0  * V0;    % mPh
         0;             % mB  (kumulierte Base startet bei 0)
         DOT0;];          % DOT

end


function s = messwert(feld, ab, name)
% zieht Zeit + Wert raus und rechnet die Messvarianz (a*y+b)^2
t = feld.BatchAge(:);
y = feld.Wert(:);
s.t   = t;
s.y   = y;
% if name == "Base"
%    s.var = feld.Variance(:)/1000; % Umrechnung von mL in L
% elseif name == "O2"
%     s.var = feld.Variance(:) * 100;
% else
% s.var = feld.Variance(:);
% end
s.var = (ab(1).*y + ab(2)).^2;
end