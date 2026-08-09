function [nRemoved, fig] = plot_cleaning(datenordner, name)
% PLOT_CLEANING  Zeigt fuer einen Datensatz, welche Punkte entfernt wurden.
%
%   [nRemoved, fig] = plot_cleaning(datenordner, name)
%
%   blau  o  = behalten
%   rot   x  = entfernt
%
% Verwendet dieselbe Funktion clean_mess.m wie das Preprocessing, damit der
% Plot garantiert das zeigt, was der Fit tatsaechlich benutzt. Das Speichern
% uebernimmt der Aufrufer (create_cleaned_plots.m).

Mess_raw = load(fullfile(datenordner, sprintf('Mess_RamScDef%s.mat', name))).Mess;
Mess_cln = clean_mess(Mess_raw, name, false);      % ohne Konsolenausgabe

% {Feldname im Mess-Struct, Titel, y-Achse}
kanaele = { 'Biomasse', 'Biomasse', 'c_X (g/L)';
            'Glucose',  'Glucose',  'c_{Glc} (g/L)';
            'Ammonium', 'Ammonium', 'c_{Am} (g/L)';
            'Phosphat', 'Phosphat', 'c_{Ph} (g/L)';
            'BASE',     'Base',     'm_B (L)';
            'O2',       'DOT',      'DOT (%)' };

fig = figure('Name', sprintf('Bereinigung RamScDef%s', name), ...
             'Position', [150 40 950 1000]);

nRemoved   = 0;
legendDone = false;

for i = 1:size(kanaele,1)
    t_raw = Mess_raw.Messdaten.(kanaele{i,1}).BatchAge(:);
    y_raw = Mess_raw.Messdaten.(kanaele{i,1}).Wert(:);
    y_cln = Mess_cln.Messdaten.(kanaele{i,1}).Wert(:);

    % Felder koennen unterschiedlich lang sein -> gemeinsamer Bereich
    n = min([numel(t_raw) numel(y_raw) numel(y_cln)]);
    t_raw = t_raw(1:n);  y_raw = y_raw(1:n);  y_cln = y_cln(1:n);

    % Entfernt = war vorher endlich, ist jetzt NaN
    removed  = isfinite(y_raw) & ~isfinite(y_cln);
    kept     = isfinite(y_cln) & isfinite(t_raw);
    nRemoved = nRemoved + nnz(removed);

    subplot(6,1,i);
    plot(t_raw(kept), y_raw(kept), 'o', 'MarkerSize', 4, ...
         'MarkerFaceColor', [0.20 0.40 0.90], ...
         'MarkerEdgeColor', [0.20 0.40 0.90]); hold on;

    if any(removed)
        plot(t_raw(removed), y_raw(removed), 'x', 'MarkerSize', 9, ...
             'LineWidth', 2, 'Color', [0.85 0.15 0.15]);
        title(sprintf('%s  --  %d entfernt', kanaele{i,2}, nnz(removed)));
        if ~legendDone
            legend('behalten','entfernt','Location','best');
            legendDone = true;
        end
    else
        title(kanaele{i,2});
    end

    ylabel(kanaele{i,3}); grid on;
end

xlabel('BatchAge (h)');
sgtitle(sprintf('RamScDef%s -- Datenbereinigung (%d Punkte entfernt)', ...
        name, nRemoved));
end