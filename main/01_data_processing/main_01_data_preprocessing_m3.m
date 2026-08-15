%% Preprocessing Modell 3 (Fed-Batch mit Ethanol)
% Training:    RamScDef10
% Validierung: RamScDef07
clear; clc; close all;


filepath = pwd;
datenordner = fullfile(filepath, 'Daten', 'MessDaten_SPI1_Projekt');

train = load(fullfile(datenordner, 'Mess_RamScDef07.mat')).Mess;   % zum Trainieren
val   = load(fullfile(datenordner, 'Mess_RamScDef10.mat')).Mess;   % zum Validieren
ab = []; 
% a/b aus dem Krämer und King 2017 Paper
% ab.Biomasse = [0.02  0.015];
% ab.Glucose  = [0.06  0.25 ];
% ab.Ammonium = [0.06  0.01 ];
% ab.Phosphat = [0.07  0.01];
% ab.O2       = [0    0.5]; % Muss nochmal recherchiert werden. Steht nicht in Krämer und King
% ab.Base     = [0.01  0.01]; % Base steht in der Doku mit 0.01 und 10 aber Einheit ist in mL angegeben und im Modell wird nur mit L gerechnet



% Training und Validierung mit derselben Routine aufbereiten
TrainData = aufbereiten(train, ab);
TrainDataProbe = train.Probenahmen;
ValData   = aufbereiten(val,   ab);
ValDataProbe = val.Probenahmen;

% Pruefen ob die Varianzen stimmen
check_var(TrainData, 'traindata')
check_var(ValData, 'valdata')


fprintf('Summe Probenvolumen = %.4g (erwartet ~0.2-0.7 L)\n', sum(TrainDataProbe.Volumen));
fprintf('max. Einzelprobe    = %.4g (erwartet ~0.01-0.03 L)\n', max(TrainDataProbe.Volumen));



% speichern
zielordner = fullfile(filepath, 'Daten', 'Daten_Processed');
if ~exist(zielordner, 'dir'); mkdir(zielordner); end
ziel = fullfile(zielordner, 'Processed_FedBatch_Modell3.mat');
save(ziel, 'TrainData', 'ValData', "TrainDataProbe", "ValDataProbe");

fprintf('Fertig. Gespeichert unter:\n%s\n', ziel);


%% ---------------------------------------------------------------
function D = aufbereiten(Mess, ab)

M = Mess.Messdaten;
% Messgroessen einlesen (Zeit in h, Wert in g/L bzw. % bei O2)
% D.Biomasse = messwert(M.Biomasse, ab.Biomasse, "Biomasse");
% D.Glucose  = messwert(M.Glucose,  ab.Glucose, "Glucose");
% D.Ammonium = messwert(M.Ammonium, ab.Ammonium, "Ammonium");
% D.Phosphat = messwert(M.Phosphat, ab.Phosphat, "Phosphat");
% D.O2       = messwert(M.O2,       ab.O2, "O2");     % = DOT in %
% D.Base     = messwert(M.BASE,     ab.Base, "Base");

% Berechnung der Varianz aus Parametern aus Messwerten und nicht aus
% Literatur
D.Biomasse = messwert(M.Biomasse, [M.Biomasse.VarParam.a M.Biomasse.VarParam.b], "Biomasse");
D.Glucose  = messwert(M.Glucose,  [M.Glucose.VarParam.a M.Glucose.VarParam.b], "Glucose");
D.Ammonium = messwert(M.Ammonium, [M.Ammonium.VarParam.a M.Ammonium.VarParam.b], "Ammonium");
D.Phosphat = messwert(M.Phosphat, [M.Phosphat.VarParam.a M.Phosphat.VarParam.b], "Phosphat");
D.O2       = messwert(M.O2,       [M.O2.VarParam.a M.O2.VarParam.b], "O2");     % = DOT in %
D.Base     = messwert(M.BASE,     [M.BASE.VarParam.a M.BASE.VarParam.b], "Base");

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
mB0 = D.Base.y(1);
DOT0  = D.O2.y(1);

D.x0 = [ V0;            % V
         cX0   * V0;    % mX
         cGlc0 * V0;    % mGlc
         cAm0  * V0;    % mAm
         cPh0  * V0;    % mPh
         mB0;             % mB  (kumulierte Base startet nicht ueberall bei 0)
         DOT0;];          % DOT

end

function s = messwert(feld, ab, name)
% zieht Zeit + Wert raus und rechnet die Messvarianz (a*y+b)^2
[t, ord] = sort(feld.BatchAge(:));
y = feld.Wert(:);  
y = y(ord); % messwerte ordnen, ein paar falsch sortierte messwerte haben probleme gemacht im time weighting

%t = feld.BatchAge(:);
%y = feld.Wert(:);

s.t   = t;
s.y   = y;

switch name
    case "Base"
        % sigma^2_mL = a*(1000*y_L) + b   [mL^2]   -> /1e6 fuer L^2
        s.var = ab(1).*y./1000 + ab(2)./1e6; % wir haben davor durch 1000 geteilt aber bei der varianz muss quadriert werden
    case "O2"
        % sigma^2_frac = a*(y_%/100) + b           -> *1e4 fuer %^2
        s.var = 100.*ab(1).*y + 1e4.*ab(2); % O2 war um einen faktor von 10 falsch
    otherwise
        s.var = ab(1).*y + ab(2);
end
s.var = max(s.var, eps);
s.var_file = feld.Variance(:);   % Varianz aus dem Datensatz
end

%if name == "Base"
%    s.var = s.var/1000; % von mL in L
%end
%if name == "O2"
%    s.var = s.var(:)*100; % Umrechnung in %
%end

%end

% pruefen ob die Varainzberechnungen stimmen
function check_var(D, label)
fprintf('--- Varianz-Check %s ---\n', label);
flds = {'Biomasse','Glucose','Ammonium','Phosphat','Base','O2'};
for i = 1:numel(flds)
    m = D.(flds{i});
    sig_ab = median(sqrt(m.var));

    % (a) empirisch aus den Daten -- nur bei dichter Abtastung sinnvoll
    if numel(m.y) >= 50
        sig_emp = std(diff(m.y,2)) / sqrt(6);
        fprintf('  %-9s: sigma_ab=%.4g  sigma_emp=%.4g  Faktor=%.2f\n', ...
                flds{i}, sig_ab, sig_emp, sig_ab/sig_emp);
    end

    % (b) gegen die mitgelieferte Varianz
    if isfield(m,'var_file') && ~isempty(m.var_file)
        sig_file = median(sqrt(m.var_file(:)), 'omitnan');
        fprintf('  %-9s: sigma_ab=%.4g  sigma_file=%.4g  Faktor=%.2f\n', ...
                flds{i}, sig_ab, sig_file, sig_ab/sig_file);
    end
end
end
