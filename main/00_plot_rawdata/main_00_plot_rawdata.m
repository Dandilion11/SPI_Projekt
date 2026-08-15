%% main_00_plot_rawdata.m
%  Stellt die ROHDATEN aller derzeit genutzten Experimente dar -
%  eine Figur je Experiment, ungeschnitten direkt aus den Mess-Dateien.
%  Relevante FEED-ZEITPUNKTE (Glucose, Base, Saeure) werden markiert;
%  pro Experiment wird eine u-Zeilen-Diagnose in die Konsole geschrieben.
%
%  Genutzte Experimente:
%    Modell 1 & 2 (Batch)    : 03, 04, 06
%    Modell 3 (Fed-Batch)    : 07, 10, 06
%  -> Vereinigung: 03, 04, 06, 07, 10
%
%  u-Konvention (Feedraten in L/h): 2=uAm, 4=uPh, 6=uGlc, 9=uAcid, 10=uBase.
%  (Zeilen 3/5/7 sind konstante Eingangs-Konzentrationen, keine Feeds;
%   Zeile 8 = Antischaum, fuer die Modelle irrelevant.)
% =========================================================================
clear; clc; close all;

scriptDir   = pwd;
datenordner = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt');

expNums = [3 4 6 7 10];

% Nur Feed-Laeufe zaehlen, deren Spitzenrate >= feedRelThresh * Max-Rate ist.
% Filtert winzige/verrauschte Nicht-Null-Eintraege heraus (z.B. Priming).
feedRelThresh = 0.05;   % 5 % der maximalen Feedrate

% Feed-Gruppen: {Anzeigename, u-Zeile(n), Farbe} - nur modellrelevante Feeds
groups = { 'Glucose',  6,      [0.85 0.33 0.10]; ...
           'Base',     10,     [0.40 0.40 0.40]; ...
           'Saeure',   9,      [0.30 0.55 0.75] };

% Kanaele: {Feldname in Messdaten, Anzeigename, y-Label}
kanal = { 'Biomasse', 'Biomasse',         'c_X (g/L)';     ...
          'Glucose',  'Glucose',          'c_{Glc} (g/L)'; ...
          'Ammonium', 'Ammonium',         'c_{Am} (g/L)';  ...
          'Phosphat', 'Phosphat',         'c_{Ph} (g/L)';  ...
          'BASE',     'Base',             'Base';          ...
          'O2',       'Sauerstoff (DOT)', 'DOT (%)'        };

for e = 1:numel(expNums)
    num   = expNums(e);
    fname = fullfile(datenordner, sprintf('Mess_RamScDef%02d.mat', num));
    if ~isfile(fname), warning('Datei fehlt: %s', fname); continue; end

    Mess = load(fname).Mess;
    M    = Mess.Messdaten;
    if isfield(Mess, 'Name'), intName = Mess.Name; else, intName = ''; end
    u = [];
    if isfield(Mess, 'u'), u = Mess.u; end

    % --- Feed-Onsets je Gruppe (schwellwertgefiltert) + Diagnose ---
    grpOnsets = cell(size(groups,1),1);
    fprintf('\n=== RamScDef%02d (intern: %s): Feed-Diagnose (>= %.0f%% Max-Rate) ===\n', ...
            num, intName, feedRelThresh*100);
    if isempty(u)
        fprintf('  keine u-Matrix vorhanden.\n');
    else
        for g = 1:size(groups,1)
            [ons, pk, nRaw] = feed_onsets(u, groups{g,2}, feedRelThresh);
            grpOnsets{g} = ons;
            if nRaw == 0
                fprintf('  %-8s (u%2d): inaktiv\n', groups{g,1}, groups{g,2});
            elseif isempty(ons)
                fprintf('  %-8s (u%2d): %d Roh-Ereignisse, alle < Schwelle (vernachlaessigbar)\n', ...
                        groups{g,1}, groups{g,2}, nRaw);
            else
                fprintf('  %-8s (u%2d): %d von %d Ereignissen relevant; erstes %.2f h; Peak %.4g L/h\n', ...
                        groups{g,1}, groups{g,2}, numel(ons), nRaw, ons(1), max(pk));
            end
        end
    end

    % --- Figur ---
    figure('Name', sprintf('Rohdaten RamScDef%02d', num), 'Position', [100 60 1000 760]);
    sgtitle(sprintf('Rohdaten RamScDef%02d   (intern: %s)', num, intName), 'Interpreter', 'none');

    for k = 1:size(kanal,1)
        subplot(3, 2, k);
        f = kanal{k,1};
        if ~isfield(M, f)
            title(sprintf('%s (nicht vorhanden)', kanal{k,2})); axis off; continue;
        end
        t = M.(f).BatchAge(:);
        y = M.(f).Wert(:);
        if isempty(y) || ~isnumeric(y) || ~isnumeric(t)
            title(sprintf('%s (keine numerischen Daten)', kanal{k,2})); axis off; continue;
        end

        dense  = numel(y) > 200;
        hasVar = isfield(M.(f), 'Variance') && numel(M.(f).Variance) == numel(y);
        if dense
            plot(t, y, '-', 'LineWidth', 1);
        elseif hasVar
            errorbar(t, y, sqrt(max(M.(f).Variance(:), 0)), 'o-', 'MarkerSize', 3);
        else
            plot(t, y, 'o-', 'MarkerSize', 3);
        end
        hold on;

        for g = 1:size(groups,1)
            for j = 1:numel(grpOnsets{g})
                xline(grpOnsets{g}(j), '--', 'Color', groups{g,3}, 'HandleVisibility', 'off');
            end
        end

        if k == 1
            hLeg = []; nam = {};
            for g = 1:size(groups,1)
                if ~isempty(grpOnsets{g})
                    hLeg(end+1) = plot(nan, nan, '--', 'Color', groups{g,3}, 'LineWidth', 1.5); %#ok<AGROW>
                    nam{end+1}  = groups{g,1}; %#ok<AGROW>
                end
            end
            if ~isempty(hLeg), legend(hLeg, nam, 'Location', 'best'); end
        end

        title(kanal{k,2}); xlabel('BatchAge (h)'); ylabel(kanal{k,3}); grid on;
    end
end

fprintf('\nRohdaten-Plots erstellt (Feeds >= %.0f%% der Max-Rate markiert).\n', feedRelThresh*100);


%% ========================================================================
function [tf, pk, nRaw] = feed_onsets(u, rows, relThresh)
% Beginn jedes zusammenhaengenden Fuetterungs-Laufs, gefiltert nach
% Spitzenrate. Ein Lauf zaehlt nur, wenn seine maximale Rate
% >= relThresh * (groesste Rate im Signal) ist.
% Rueckgabe: gefilterte Onset-Zeiten tf, deren Spitzenraten pk,
% und die Zahl aller (ungefilterten) Roh-Laeufe nRaw.
    tf = []; pk = []; nRaw = 0;
    rows = rows(rows <= size(u,1));
    if isempty(rows), return; end

    tu  = u(1,:);
    sig = max(u(rows,:), [], 1);      % groesste Rate ueber die Gruppe
    act = sig > 0;
    d   = diff([0 act 0]);
    s   = find(d == 1);               % Lauf-Starts
    en  = find(d == -1) - 1;          % Lauf-Enden
    nRaw = numel(s);
    if nRaw == 0, return; end

    gmax = max(sig);
    for i = 1:nRaw
        p = max(sig(s(i):en(i)));
        if gmax <= 0 || p >= relThresh * gmax
            tf(end+1) = tu(s(i)); %#ok<AGROW>
            pk(end+1) = p;        %#ok<AGROW>
        end
    end
end
