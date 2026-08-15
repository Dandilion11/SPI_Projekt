% main_01_data_preprocessing_m3_extended.m
%
% Preprocessing Modell 3 (Fed-Batch, ohne Ethanol) ueber mehrere Experimente.
%   Training:    RamScDef03, 04, 06, 10
%   Validierung: RamScDef07
%   Def05 faellt raus (keine O2-Messung).
%
% Warum diese Aufteilung:
%   Def04 bleibt durchgehend ueber 60 % DOT, also ohne Sauerstofflimitierung
%   -- der Bereich, in dem die Modellklasse gilt. Def06 hat 15 L Startvolumen
%   statt 10 L und laeuft bei t~30 h in Ammoniummangel. Def10 liefert den
%   laengsten Horizont und die hoechste Biomasse. Def07 bleibt als
%   Validierungslauf zurueckgehalten.
%
% rpm ist in allen Laeufen 600, deshalb ist ein gemeinsames KLa zulaessig.
%
% Die Bereinigung liegt in clean_mess.m, damit dieses Skript und
% plot_cleaning.m garantiert dieselbe Logik verwenden.

clear; clc; close all;

filepath    = pwd;
datenordner = fullfile(filepath, 'Daten','MessDaten_SPI1_Projekt');

trainNames = {'03','04','06','10'};
valNames   = {'07'};

%% Laden, bereinigen, aufbereiten ----------------------------------------
fprintf('=== DATENBEREINIGUNG ===\n');

TrainSet = struct([]);
for k = 1:numel(trainNames)
    TrainSet = append_experiment(TrainSet, datenordner, trainNames{k});
end
ValSet = struct([]);
for k = 1:numel(valNames)
    ValSet = append_experiment(ValSet, datenordner, valNames{k});
end

%% Kontrolle: hat die Bereinigung gegriffen? -----------------------------
% Wenn BatchAge und Wert unterschiedlich lang sind, kann eine Zuweisung ins
% Leere laufen, ohne dass MATLAB warnt
i03 = find(strcmp({TrainSet.name},'RamScDef03'), 1);
fprintf('  Def03 DOT_max = %.1f %%  (Soll <= 100)\n', max(TrainSet(i03).O2.y));

%% Uebersicht und Varianz-Check ------------------------------------------
fprintf('\n=== UEBERSICHT TRAINING ===\n');    report_set(TrainSet);
fprintf('\n=== UEBERSICHT VALIDIERUNG ===\n'); report_set(ValSet);

fprintf('\n=== VARIANZ-CHECK ===\n');
for k = 1:numel(TrainSet), check_var(TrainSet(k)); end
for k = 1:numel(ValSet),   check_var(ValSet(k));   end

%% Speichern --------------------------------------------------------------
zielordner = fullfile(filepath,'Daten','Daten_Processed');
if ~exist(zielordner,'dir'); mkdir(zielordner); end
ziel = fullfile(zielordner,'Processed_FedBatch_Modell3_MultiExp.mat');
save(ziel, 'TrainSet', 'ValSet');

fprintf('\nGespeichert: %s\n', ziel);


%% ======================================================================
%  Hilfsfunktionen
%% ======================================================================

function S = append_experiment(S, datenordner, name)
% Laedt einen Datensatz, reinigt ihn und haengt ihn ans Struct-Array an.
    f = fullfile(datenordner, sprintf('Mess_RamScDef%s.mat', name));
    Mess = load(f).Mess;
    Mess = clean_mess(Mess, name);
    D = aufbereiten(Mess);
    D.name = ['RamScDef' name];
    if isempty(S), S = D; else, S(end+1) = D; end
end


function D = aufbereiten(Mess)
% Baut aus einem Mess-Struct die im Fit verwendete Datenstruktur.
    M = Mess.Messdaten;

    % Messkanal im Ergebnis -> Feldname im Mess-Struct
    kanaele = {'Biomasse','Biomasse'; 'Glucose','Glucose'; ...
               'Ammonium','Ammonium'; 'Phosphat','Phosphat'; ...
               'O2','O2';             'Base','BASE'};
    for i = 1:size(kanaele,1)
        D.(kanaele{i,1}) = messwert(M.(kanaele{i,2}), kanaele{i,1});
    end

    % Volumen (Dichte 1 kg/L angenommen -> V[L] = Gewicht[kg]).
    % In Def03/04/06 ist Weight nur ein Skalar, nur Def07/10 haben eine
    % echte Zeitreihe.
    vt = M.Weight.BatchAge(:);  vy = M.Weight.Wert(:);
    n  = min(numel(vt), numel(vy));
    ok = isfinite(vt(1:n)) & isfinite(vy(1:n));
    D.V.t = vt(ok);  D.V.y = vy(ok);
    D.V.isTimeSeries = nnz(ok) > 5;
    V0 = D.V.y(1);

    D.u     = Mess.u;              % Stellgroessen
    D.Probe = Mess.Probenahmen;    % Probenahmezeiten und -volumina

    % Startwerte, Reihenfolge wie in Modell3_woEtOH_10p:
    % x = [V; mX; mGlc; mAm; mPh; mB; DOT]. Die Messungen sind
    % Konzentrationen, das Modell rechnet in Massen -> mal V0.
    D.x0 = [ V0;
             D.Biomasse.y(1) * V0;
             D.Glucose.y(1)  * V0;
             D.Ammonium.y(1) * V0;
             D.Phosphat.y(1) * V0;
             D.Base.y(1);                  % kumuliert, startet nicht bei 0
             D.O2.y(1) ];

    % x0 gilt formal bei t = u(1,1). Liegt die erste Messung deutlich
    % spaeter, ist die Anfangsbedingung inkonsistent -> report_set warnt.
    tfs = [D.Biomasse.t(1) D.Glucose.t(1) D.Ammonium.t(1) ...
           D.Phosphat.t(1) D.Base.t(1)     D.O2.t(1)];
    D.x0_offset = max(abs(tfs - D.u(1,1)));
end


function s = messwert(feld, name)
% Zeit und Wert extrahieren, NaN entfernen, sortieren, Varianz berechnen.
%
% Sortieren ist zwingend: x0 nimmt y(1), also den ersten Wert NACH der
% Zeit -- nicht den ersten im Speicher.
%
% Varianzmodell sigma^2 = a*y + b mit a,b aus dem Datensatz. a und b sind
% in den Doku-Einheiten kalibriert, y liegt hier schon in Modell-Einheiten
% vor. Varianz skaliert QUADRATISCH mit einer Einheitenumrechnung:
%   Base: Doku mL, Modell L     -> a*y/1000 + b/1e6
%   O2:   Doku Anteil, Modell % -> 100*a*y  + 1e4*b

    t = feld.BatchAge(:);
    y = feld.Wert(:);
    vf = nan(size(y));
    if isfield(feld,'Variance') && ~isempty(feld.Variance)
        m = min(numel(vf), numel(feld.Variance));
        vf(1:m) = feld.Variance(1:m);
    end

    % Felder koennen unterschiedlich lang sein -> gemeinsamer Bereich
    n = min([numel(t) numel(y) numel(vf)]);
    t = t(1:n);  y = y(1:n);  vf = vf(1:n);

    ok = isfinite(t) & isfinite(y);          % u.a. die NaN aus clean_mess
    t = t(ok);  y = y(ok);  vf = vf(ok);
    [t, ord] = sort(t);  y = y(ord);  vf = vf(ord);

    s.t = t;
    s.y = y;

    a = feld.VarParam.a;  b = feld.VarParam.b;
    switch name
        case 'Base', s.var = a.*y/1000 + b/1e6;
        %case 'O2',   s.var = 100*a.*y  + 1e4*b;
        otherwise,   s.var = a.*y + b;
    end
    s.var      = max(s.var, eps);
    s.var_file = vf;      % unkonvertierte Varianz, nur fuer check_var
end


function report_set(S)
% Kurzuebersicht je Experiment: Umfang, Startwerte, Feeds, Probenahmen.
    for k = 1:numel(S)
        D = S(k);  u = D.u;
        idxFeed = find(any(u([2 4 6],:) > 0, 1), 1, 'first');
        if isempty(idxFeed), t_feed = NaN; else, t_feed = u(1, idxFeed); end

        fprintf('--- %s ---\n', D.name);
        fprintf('  t: %.2f .. %.2f h | V0=%.2f L (Zeitreihe: %d) | erster Feed t=%.2f h\n', ...
            min(D.Biomasse.t), max(D.Biomasse.t), D.x0(1), D.V.isTimeSeries, t_feed);
        fprintf('  cX: %.2f .. %.2f g/L | cGlc_max=%.1f | cAm_min=%.3f | DOT: %.1f .. %.1f %%\n', ...
            min(D.Biomasse.y), max(D.Biomasse.y), max(D.Glucose.y), ...
            min(D.Ammonium.y), min(D.O2.y), max(D.O2.y));
        fprintf('  n: Bio=%d Glc=%d Am=%d Ph=%d Base=%d O2=%d\n', ...
            numel(D.Biomasse.t), numel(D.Glucose.t), numel(D.Ammonium.t), ...
            numel(D.Phosphat.t), numel(D.Base.t), numel(D.O2.t));
        fprintf('  Probenahmen: n=%d, Summe=%.3f L\n', ...
            numel(D.Probe.Volumen), sum(D.Probe.Volumen));
        if D.x0_offset > 0.5
            fprintf('  [WARNUNG] erste Messung bis zu %.2f h nach u(1,1) -> x0 inkonsistent\n', ...
                D.x0_offset);
        end
    end
end


function check_var(D)
% Plausibilitaet der Varianz: berechnetes sigma gegen (a) das empirische
% Rauschen aus zweiten Differenzen und (b) die mitgelieferte Varianz.
% Faktoren nahe 1 heissen: das Varianzmodell passt zu den Daten.
    fprintf('--- Varianz-Check %s ---\n', D.name);
    for f = {'Biomasse','Glucose','Ammonium','Phosphat','Base','O2'}
        m = D.(f{1});
        if isempty(m.y), continue; end
        sig_ab = median(sqrt(m.var));

        % (a) nur bei dichter Abtastung aussagekraeftig, sonst dominiert
        %     die echte Prozesskruemmung das Ergebnis
        if numel(m.y) >= 50
            sig_emp = std(diff(m.y,2)) / sqrt(6);
            fprintf('  %-9s sigma_ab=%.4g  sigma_emp=%.4g   Faktor=%.2f\n', ...
                    f{1}, sig_ab, sig_emp, sig_ab/sig_emp);
        end
        % (b) gegen die unkonvertierte Varianz aus dem Datensatz
        if any(isfinite(m.var_file))
            sig_file = median(sqrt(m.var_file(isfinite(m.var_file))));
            fprintf('  %-9s sigma_ab=%.4g  sigma_file=%.4g  Faktor=%.2f\n', ...
                    f{1}, sig_ab, sig_file, sig_ab/sig_file);
        end
    end
end