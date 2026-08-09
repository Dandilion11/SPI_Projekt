% create_cleaned_plots.m
%
% Erzeugt fuer jeden verwendeten Datensatz eine Abbildung, die zeigt, welche
% Messpunkte in der Bereinigung entfernt wurden.
%   blau  o  = behalten
%   rot   x  = entfernt
%
% Benoetigt clean_mess.m (die Bereinigungslogik, wird auch vom Preprocessing
% verwendet) und plot_cleaning.m (die Abbildung).

clear; clc; close all;

filepath    = pwd;
datenordner = fullfile(filepath,'Daten','MessDaten_SPI1_Projekt');
bildordner  = fullfile(filepath,'Bilder','Bereinigung');

alleNamen = {'03','04','06','07','10'};   % 05 raus (keine O2-Messung)

fprintf('=== BEREINIGUNGS-ABBILDUNGEN ===\n');

gesamt = 0;
for k = 1:numel(alleNamen)
    name = alleNamen{k};
    fprintf('RamScDef%s ...\n', name);
    [n, fig] = plot_cleaning(datenordner, name);
    gesamt = gesamt + n;
    fprintf('  %d Punkte entfernt\n', n);
    save_fig(fig, sprintf('Bereinigung_RamScDef%s', name), bildordner);
end

fprintf('\nFertig. Insgesamt %d Punkte entfernt ueber %d Datensaetze.\n', ...
        gesamt, numel(alleNamen));
fprintf('Abbildungen in: %s\n', bildordner);


%% ======================================================================
function save_fig(fig, name, ordner)
% Speichert eine Figure hell (fuer Folien) als PNG und PDF.
    if ~exist(ordner,'dir'), mkdir(ordner); end
    ziel = fullfile(ordner, [name '.png']);

    try                                   % ab R2025a: helles Theme
        set(fig, 'Theme', 'light');  drawnow;
    catch                                 % sonst von Hand umfaerben
        set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
        for a = findall(fig,'Type','axes').'
            set(a, 'Color','w', 'XColor','k', 'YColor','k', ...
                   'GridColor',[0.15 0.15 0.15], 'GridAlpha',0.15);
            set([a.Title a.XLabel a.YLabel], 'Color', 'k');
        end
        set(findall(fig,'Type','legend'), 'TextColor','k', 'Color','w');
        set(findall(fig,'Type','text'), 'Color', 'k');
        drawnow;
    end

    exportgraphics(fig, ziel, 'Resolution', 300);
    exportgraphics(fig, strrep(ziel,'.png','.pdf'), 'ContentType','vector');
    fprintf('  gespeichert: %s\n', ziel);
end