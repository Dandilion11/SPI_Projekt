clear;    
clc;       
close all;  

%% 2. Pfad-Management und Daten laden
scriptDir = fileparts(mfilename('fullpath'));
ordnerName = fullfile(scriptDir, '..', 'Daten', 'MessDaten_SPI1_Projekt'); 
dateien = dir(fullfile(ordnerName, '*.mat'));
alleExperimente = cell(1, length(dateien));


for i = 1:length(dateien)
    aktuellerPfad = fullfile(ordnerName, dateien(i).name);
    geladeneDaten = load(aktuellerPfad);

    % Daten laden, Speichern des Inhalts von 'Mess' in die jeweilige Zelle
    alleExperimente{i} = geladeneDaten.Mess; 

    fprintf('Datei %s erfolgreich geladen.\n', dateien(i).name);
    messdaten = load(aktuellerPfad);

end

%% 3. Daten plotten – ein Figure pro Experiment
for i = 1:length(alleExperimente)
    exp_i  = alleExperimente{i};
    felder = fieldnames(exp_i.Messdaten);

    nFelder = length(felder);
    nCols   = 3;
    nRows   = ceil((nFelder + 1) / nCols);

    figure('Name', exp_i.Name, 'NumberTitle', 'off');
    sgtitle(exp_i.Name, 'Interpreter', 'none');

    % --- Messdaten-Felder ---
    for k = 1:nFelder
        feld = felder{k};
        t_k  = exp_i.Messdaten.(feld).BatchAge;
        y_k  = exp_i.Messdaten.(feld).Wert;

        subplot(nRows, nCols, k);

        if isscalar(t_k) || isscalar(y_k)
            yline(double(y_k), 'r', 'LineWidth', 1.5);
        else
            stem(t_k(:), y_k(:), 'filled', 'MarkerSize', 4);
        end

        title(feld, 'Interpreter', 'none');
        xlabel('BatchAge (h)');
        grid on;

        if isfield(exp_i.Messdaten.(feld), 'Variance')
            var_k = exp_i.Messdaten.(feld).Variance;
            if numel(var_k) == numel(y_k) && numel(y_k) > 1
                hold on;
                errorbar(t_k(:), y_k(:), sqrt(abs(var_k(:))), ...
                    'LineStyle', 'none', 'Color', [0.7 0.7 0.7]);
                hold off;
            end
        end
    end

    % --- Eingang u gegen t (robust) ---
    subplot(nRows, nCols, nFelder + 1);

    t_u = exp_i.t(:);        % Spaltenvektor% Spaltenvektor
    u   = exp_i.u;

    % u so ausrichten dass erste Dimension = length(t_u)
    if size(u, 1) == length(t_u)
        % passt direkt
    elseif size(u, 2) == length(t_u)
        u = u';              % transponieren
    else
        % Dimensionen passen nicht – Fallback: nur Index
        fprintf('Experiment %d: t(%d) passt nicht zu u(%dx%d) – plotte gegen Index\n', ...
            i, length(t_u), size(exp_i.u,1), size(exp_i.u,2));
        t_u = (1:size(u,2))';
        if size(u,1) ~= length(t_u), u = u'; end
    end

    plot(t_u, u);
    title('Eingang u', 'Interpreter', 'none');
    xlabel('t (h)');
    grid on;
    legendStr = arrayfun(@(x) sprintf('u_%d',x), 1:size(u,2), ...
        'UniformOutput', false);
    legend(legendStr, 'Location', 'best');
end