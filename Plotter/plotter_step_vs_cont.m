clear;    
clc;       
close all;

%% 1. Konfiguration
scriptDir = fileparts(mfilename('fullpath'));
ordnerName = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt'); 
dateien = dir(fullfile(ordnerName, '*.mat'));


if length(dateien) < 6
    warning('Weniger als 6 Dateien gefunden. Überprüfe den Ordner %s.', ordnerName);
end

% Definierte Variablen für den Plot
vars_to_plot = {'Biomasse', 'Glucose', 'Phosphat', 'Ammonium', ...
                'BASE', 'O2', 'm_Feed_C', 'm_Feed_Am', 'm_Feed_Ph'};

% Gruppendefinition (Indizes der Dateien basierend auf der alphabetischen Sortierung)
gruppen = {1:3, 4:6};
gruppen_namen = {'Gruppe 1: Pulse-Feeding (Exp 1-3)', 'Gruppe 2: Fed-Batch (Exp 4-6)'};

% Farbmatrix für bis zu 3 Experimente pro Plot
farben = [0 0.4470 0.7410; ... % Blau
          0.8500 0.3250 0.0980; ... % Orange
          0.9290 0.6940 0.1250];    % Gelb

%% 2. Daten laden und gruppiert plotten
for g = 1:length(gruppen)
    idx_group = gruppen{g};
    
    % Überprüfen, ob die Indizes die Anzahl der Dateien nicht überschreiten
    idx_group = idx_group(idx_group <= length(dateien)); 
    if isempty(idx_group), continue; end
    
    % Neues Figure pro Gruppe erstellen
    figure('Name', gruppen_namen{g}, 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);
    sgtitle(gruppen_namen{g}, 'Interpreter', 'none', 'FontWeight', 'bold');
    
    % Array für Subplot-Achsen initialisieren (3x3 Raster für 9 Variablen)
    ax = zeros(1, length(vars_to_plot));
    for v = 1:length(vars_to_plot)
        ax(v) = subplot(3, 3, v);
        hold(ax(v), 'on'); 
        title(ax(v), vars_to_plot{v}, 'Interpreter', 'none');
        xlabel(ax(v), 'BatchAge (h)'); 
        grid(ax(v), 'on');
    end
    
    % Schleife über die Experimente in der aktuellen Gruppe
    for i_idx = 1:length(idx_group)
        i = idx_group(i_idx);
        aktuellerPfad = fullfile(ordnerName, dateien(i).name);
        geladeneDaten = load(aktuellerPfad);
        exp_i = geladeneDaten.Mess;
        
        exp_name = exp_i.Name;
        aktuelle_farbe = farben(i_idx, :);
        
        % Schleife über die auszuwählenden Variablen
        for v = 1:length(vars_to_plot)
            var_name = vars_to_plot{v};
            
            % Prüfen, ob die Variable in den Messdaten existiert
            if isfield(exp_i.Messdaten, var_name)
                t_k = exp_i.Messdaten.(var_name).BatchAge;
                y_k = exp_i.Messdaten.(var_name).Wert;
                
                % Leere Variablen überspringen (z.B. O2 in RamScDef05)
                if isempty(t_k) || isempty(y_k) || numel(y_k) == 0
                    continue; 
                end
                
                % Plot-Entscheidung: Kontinuierlich (>50 Punkte) vs. Diskret (Offline-Samples)
                is_continuous = numel(y_k) > 50;
                
                if is_continuous
                    plot(ax(v), t_k, y_k, '-', 'LineWidth', 1.5, ...
                        'Color', aktuelle_farbe, 'DisplayName', exp_name);
                else
                    % Prüfen, ob Varianz für Errorbars existiert
                    if isfield(exp_i.Messdaten.(var_name), 'Variance') && numel(exp_i.Messdaten.(var_name).Variance) == numel(y_k)
                        var_k = exp_i.Messdaten.(var_name).Variance;
                        errorbar(ax(v), t_k, y_k, sqrt(abs(var_k)), 'o--', ...
                            'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerFaceColor', aktuelle_farbe, ...
                            'Color', aktuelle_farbe, 'DisplayName', exp_name);
                    else
                        plot(ax(v), t_k, y_k, 'o--', 'LineWidth', 1.2, 'MarkerSize', 5, ...
                            'MarkerFaceColor', aktuelle_farbe, 'Color', aktuelle_farbe, 'DisplayName', exp_name);
                    end
                end
            end
        end
    end
    
    % Legende nur im ersten Subplot anzeigen (gilt für alle Variablen der Gruppe)
    legend(ax(1), 'show', 'Location', 'best');
end