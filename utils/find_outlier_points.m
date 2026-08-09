% find_outlier_points.m
%
% Durchsucht alle Rohdatensaetze systematisch nach Ausreisser-KANDIDATEN.
% Das Skript entfernt nichts -- es liefert die Belege, auf deren Grundlage
% in clean_mess.m entschieden wird.
%
% Laeuft auf den ROHEN Mess-Structs (vor aufbereiten), weil Probenahmen und
% die Stellgroessenmatrix u gebraucht werden.
%
% Grundsatz: entfernt wird nur aus modellunabhaengigen Gruenden --
% physikalisch unmoeglich, Sensorausfall, Massenbilanz-Verletzung, oder
% inkonsistent mit dem Probenahme-Protokoll. NICHT wegen grosser Residuen:
% Modell 3 ist bekanntermassen unvollstaendig (kein Ethanol), residuen-
% basiertes Putzen wuerde genau die Daten loeschen, die das belegen.
%
% Die Meldungen [ENDE] und [HINWEIS] sind ausdruecklich nur Hinweise. Beide
% schlagen auch bei echter Prozessdynamik an -- der [ENDE]-Treffer in Def10
% wurde geprueft und die Probe BEHALTEN (die Wachstumsrate ueber die 14 h
% Luecke liegt im Bereich der vorangehenden Intervalle).

clear; clc; close all;

filepath    = pwd;
datenordner = fullfile(filepath,'Daten','MessDaten_SPI1_Projekt');
alleNamen   = {'03','04','06','07','10'};   % 05 raus (keine O2-Messung)

for k = 1:numel(alleNamen)
    f = fullfile(datenordner, sprintf('Mess_RamScDef%s.mat', alleNamen{k}));
    Mess = load(f).Mess;
    scan(Mess, ['RamScDef' alleNamen{k}]);
end

fprintf('\n=== Scan abgeschlossen ===\n');


%% ======================================================================
function scan(Mess, name)
fprintf('\n========== %s ==========\n', name);
md    = Mess.Messdaten;
found = 0;
offline = {'Biomasse','Glucose','Ammonium','Phosphat'};

%% (a) Physikalisch unmoegliche Werte ------------------------------------
bad = md.O2.Wert > 100 | md.O2.Wert < 0;
if any(bad)
    fprintf('[PHYS] O2 ausserhalb [0,100]: %d Punkte (min=%.1f max=%.1f)\n', ...
            nnz(bad), min(md.O2.Wert), max(md.O2.Wert));
    found = found + 1;
end
for f = offline
    w = md.(f{1}).Wert;
    if any(w < 0)
        fprintf('[PHYS] %s negativ: %d Punkte (min=%.3f)\n', f{1}, nnz(w<0), min(w));
        found = found + 1;
    end
end

%% (b) Sensorausfall: laengere exakt konstante O2-Bloecke ----------------
[t2, o2] = sortiert(md.O2.BatchAge, md.O2.Wert);
if numel(o2) > 1
    d = [1; diff(o2)];  run0 = 0;
    for i = 1:numel(d)
        if abs(d(i)) < 1e-9
            run0 = run0 + 1;
        else
            if run0 >= 5
                fprintf('[SENSOR] O2 konstant ueber %d Punkte, t=%.2f..%.2f h (Wert %.2f)\n', ...
                        run0, t2(i-run0), t2(i-1), o2(i-1));
                found = found + 1;
            end
            run0 = 0;
        end
    end
end

%% (c) Massenbilanz: Glucose steigt ohne Feed ----------------------------
[tg, yg] = sortiert(md.Glucose.BatchAge, md.Glucose.Wert);
u = Mess.u;
for i = 2:numel(tg)
    if yg(i) - yg(i-1) > 3                       % g/L Anstieg
        win = u(1,:) >= tg(i-1) & u(1,:) <= tg(i);
        if ~any(u(6,win) > 0)
            fprintf('[BILANZ] Glucose %.1f -> %.1f g/L (t=%.2f..%.2f h) ohne uGlc>0\n', ...
                    yg(i-1), yg(i), tg(i-1), tg(i));
            found = found + 1;
        end
    end
end

%% (d) Messzeitpunkte ohne zugehoerige Probenahme ------------------------
tp = Mess.Probenahmen.BatchAge(:);
tp = tp(isfinite(tp));
for f = offline
    tb = md.(f{1}).BatchAge(:);
    tb = tb(isfinite(tb));
    for i = 1:numel(tb)
        dmin = min(abs(tp - tb(i)));
        if dmin > 0.5
            fprintf('[PROBE] %s t=%.2f h ohne Probenahme (naechste %.2f h entfernt)\n', ...
                    f{1}, tb(i), dmin);
            found = found + 1;
        end
    end
end

%% HINWEIS: letzter Messpunkt springt ------------------------------------
% Vergleicht den letzten Schritt gegen den typischen Schritt. Achtung: bei
% ungleichmaessiger Abtastung ist das kein fairer Vergleich -- ein langer
% Abstand ergibt zwangslaeufig einen grossen Schritt.
for f = offline
    [t, y] = sortiert(md.(f{1}).BatchAge, md.(f{1}).Wert);
    if numel(y) < 3, continue; end
    typ = median(abs(diff(y(1:end-1))));
    if typ > 0 && abs(y(end)-y(end-1)) > 5*typ
        fprintf('[ENDE] %s: %.2f -> %.2f bei t=%.2f h (typ. Schritt %.2f, dt=%.1f h)\n', ...
                f{1}, y(end-1), y(end), t(end), typ, t(end)-t(end-1));
        found = found + 1;
    end
end

%% HINWEIS: Einzelausreisser gegen den lokalen Median --------------------
% Robuster z-Score im gleitenden Fenster. Echte Prozessdynamik
% (Feed-Spruenge) sieht genauso aus -> niemals automatisch loeschen.
for f = offline
    [t, y] = sortiert(md.(f{1}).BatchAge, md.(f{1}).Wert);
    if numel(y) < 7, continue; end
    for i = 3:numel(y)-2
        loc = y(max(1,i-3):min(numel(y),i+3));
        loc(loc == y(i)) = [];
        if isempty(loc), continue; end
        mad_ = median(abs(loc - median(loc)));
        if mad_ > 0 && abs(y(i) - median(loc)) > 6*1.4826*mad_
            fprintf('[HINWEIS] %s t=%.2f h: %.2f (lokaler Median %.2f)\n', ...
                    f{1}, t(i), y(i), median(loc));
            found = found + 1;
        end
    end
end

if found == 0
    fprintf('  keine Auffaelligkeiten gefunden.\n');
end
end


%% ======================================================================
function [t, y] = sortiert(t, y)
% NaN entfernen und nach Zeit sortieren (Spaltenform).
    t = t(:);  y = y(:);
    n = min(numel(t), numel(y));
    t = t(1:n);  y = y(1:n);
    ok = isfinite(t) & isfinite(y);
    t = t(ok);  y = y(ok);
    [t, ord] = sort(t);  y = y(ord);
end