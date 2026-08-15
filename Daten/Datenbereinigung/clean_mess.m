function [Mess, stats] = clean_mess(Mess, name, verbose)
% CLEAN_MESS  Datensatz-spezifische Bereinigung der Messdaten.
%
%   [Mess, stats] = clean_mess(Mess, name)          % mit Konsolenausgabe
%   [Mess, stats] = clean_mess(Mess, name, false)   % ohne
%
% Eigene Datei (keine lokale Funktion im Preprocessing), damit Fit und
% Bereinigungs-Plot garantiert dieselbe Logik verwenden.
%
% GRUNDSATZ
%   Ein Messpunkt wird nur entfernt, wenn es dafuer einen Grund gibt, der
%   UNABHAENGIG vom Modell ist:
%     (a) physikalisch unmoeglicher Wert (z.B. DOT > 100 %)
%     (b) Sensorausfall (Block exakt konstant oder auf 0)
%     (c) Massenbilanz-Verletzung (Substrat erscheint ohne Feed)
%     (d) Inkonsistenz mit dem Probenahme-Protokoll
%
%   NICHT zulaessig: Entfernen wegen grosser Residuen. Modell 3 ist
%   bekanntermassen unvollstaendig (Ethanol wurde nie gemessen, deshalb
%   fehlt die Overflow-Kinetik). Residuenbasiertes Putzen wuerde genau die
%   Punkte loeschen, die diese Modellgrenze belegen.
%
% GEPRUEFT, ABER NICHTS ENTFERNT
%   (c) Zwei Glucose-Anstiege (Def03 bei t~23 h und t~48 h) sind
%       bilanzkonform -- in beiden Faellen steht ein uGlc-Feed davor.
%   (d) Alle Offline-Messzeiten liegen exakt auf Probenahmezeitpunkten
%       (max. Abweichung 0.00 h in allen fuenf Datensaetzen).
%
% BEWUSST NICHT ENTFERNT
%   - DOT-Einbrueche auf 0 (Def06/07/10): reproduzierbar, in mehreren
%     Laeufen, korreliert mit den Feed-Phasen. Ohne O2-Limitierung in rX
%     kann das Modell sie nicht abbilden -- das ist ein BEFUND, keine
%     Stoerung.
%   - Biomassewachstum bei cGlc = 0 (Def04, Def06, Def10): die Signatur
%     des Overflow-Metabolismus.
%
% Entfernte Werte werden auf NaN gesetzt, messwert() filtert sie heraus.
%
% FELDLAENGEN
%   In manchen Datensaetzen sind BatchAge, Wert und Variance unterschiedlich
%   lang. Ein Index aus BatchAge darf deshalb nicht ungeprueft auf Wert
%   angewendet werden -- MATLAB wuerde Wert stillschweigend verlaengern
%   statt zu warnen. Deshalb ueberall nur der gemeinsame Bereich.

if nargin < 3, verbose = true; end

md = Mess.Messdaten;
n0 = count_valid(md);

switch name

    case '03'
        % (a)+(b) O2-Sensorausfall im Fenster t = 18.5 .. 23.6 h.
        % Ein zusammenhaengendes Sensorereignis, kein echter DOT-Verlauf:
        %   t=18.4 h: 85.12 %   letzter gesunder Wert
        %   t=18.5 h:  2.06 %   Abfall um 83 Prozentpunkte in 6 min --
        %                       echte O2-Limitierung faellt ueber Stunden
        %                       ab (vgl. Def06/07/10), nicht in einem Schritt
        %   t=19.0 .. 22.5 h:   exakt 0.00 ueber 3.5 h
        %   t=22.8 .. 23.5 h:   98.9 -> 123.1 %. DOT ist Prozent der
        %                       Saettigung, > 100 % ist unmoeglich und
        %                       beweist den Sensorfehler.
        %   t=23.7 h: 78.66 %   wieder plausibel
        % Entfernt wird das GESAMTE Fenster, nicht nur die exakten Nullen:
        % die Punkte dazwischen gehoeren zum selben Ausfall.
        [tO, yO] = common_range(md.O2.BatchAge, md.O2.Wert);
        bad = (tO >= 18.5 & tO <= 23.6) | (yO > 100);
        md.O2.Wert(bad) = NaN;
        md.O2.Wert(1) = NaN;

    case '04'
        md.O2.Wert(1) = NaN;
        
    case {'06','07', '10'}
        % Keine Befunde nach den Kriterien (a)-(d).

    otherwise
        warning('clean_mess: unbekannter Datensatz "%s" -- keine Bereinigung.', name);
end

Mess.Messdaten = md;

%% Bilanz: wie viele Punkte wurden je Kanal entfernt? --------------------
n1    = count_valid(md);
flds  = fieldnames(n0);
stats = struct();
for i = 1:numel(flds)
    stats.(flds{i}) = n0.(flds{i}) - n1.(flds{i});
end

if verbose
    entfernt = struct2array(stats);
    if any(entfernt > 0)
        fprintf('  RamScDef%s:\n', name);
        for i = 1:numel(flds)
            if stats.(flds{i}) > 0
                fprintf('    %-9s %4d -> %4d  (%d entfernt)\n', flds{i}, ...
                        n0.(flds{i}), n1.(flds{i}), stats.(flds{i}));
            end
        end
    else
        fprintf('  RamScDef%s: keine Punkte entfernt\n', name);
    end
end
end


%% ======================================================================
function [t, y] = common_range(t, y)
% Kuerzt beide Vektoren auf die gemeinsame Laenge (Spaltenform).
    t = t(:);  y = y(:);
    n = min(numel(t), numel(y));
    t = t(1:n);  y = y(1:n);
end


function n = count_valid(md)
% Zaehlt die endlichen Messwerte je Kanal (fuer die Bilanz oben).
    n.Biomasse = nnz(isfinite(md.Biomasse.Wert));
    n.Glucose  = nnz(isfinite(md.Glucose.Wert));
    n.Ammonium = nnz(isfinite(md.Ammonium.Wert));
    n.Phosphat = nnz(isfinite(md.Phosphat.Wert));
    n.Base     = nnz(isfinite(md.BASE.Wert));
    n.O2       = nnz(isfinite(md.O2.Wert));
end