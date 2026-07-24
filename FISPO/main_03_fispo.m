%Executes the STRIKE-GOLDD toolbox for structural observability and identifiability (Task 3)

%% main_03_fispo

clear; clc; close all;


%% 1. STRIKE-GOLDD Toolbox zum Pfad hinzufügen

% Pfad zur Toolbox
scriptDir = pwd;

% Einen Ordner nach oben gehen
sg_path = fullfile(fileparts(pwd), 'STRIKE-GOLDD');

% Pfad zum Modell-Ordner
model_path = fullfile(fileparts(pwd), 'Modelle');


% STRIKE-GOLDD zum MATLAB-Pfad hinzufügen
if exist(sg_path, 'dir')
    addpath(genpath(sg_path));   % genpath fügt auch Unterordner hinzu
    fprintf('[OK] STRIKE-GOLDD gefunden:\n     %s\n', sg_path);
else
    error(['[FEHLER] STRIKE-GOLDD nicht gefunden unter:\n  %s\n\n', ...
           'Erwartete Struktur: SPI_Projekt/STRIKE-GOLDD/strike_goldd.m\n'], sg_path);
end

% Modell-Ordner zum MATLAB-Pfad hinzufügen
if exist(model_path, 'dir')
    addpath(model_path);
    fprintf('[OK] Modell-Ordner gefunden:\n     %s\n', model_path);
else
    error('[FEHLER] Modell-Ordner nicht gefunden unter:\n  %s\n\n', model_path);
end



%% 2. FISPO-Analyse: Modell 1 (Batch, Monod, ohne O2)

fprintf('\n');
fprintf('=============================================================\n');
fprintf('  FISPO-Analyse: MODELL 1 — Batch, Monod\n');
fprintf('  Zustände  : cX, cGlc\n');
fprintf('  Parameter : mumax, KS, YXS\n');
fprintf('=============================================================\n\n');

% Modell aufstellen
syms cX cGlc real       
syms mumax KS YXS real  

x        = [cX; cGlc];        % 1. Zustandsvektor
p        = [mumax; KS; YXS];  % 2. Parametervektor 
u        = [];                % 3. Bekannte Eingänge 
w        = [];                % 4. Unbekannte Eingänge / Störgrößen
ic       = [];                % 5. Anfangsbedingungen
known_ic = [];                % 6. Bekannte ICs

% Wachstumsrate und Differentialgleichungen (f)
mu = mumax * cGlc / (KS + cGlc);
f  = [ mu * cX;               
      -1/YXS * mu * cX];

% Messgleichungen (h)
h  = [cX; cGlc];              

% Als .mat speichern für Toolbox
mat_filename = fullfile(model_path, 'SF_Modell1.mat');
save(mat_filename, 'x', 'f', 'h', 'u', 'p', 'w', 'ic', 'known_ic');
fprintf('Modell-Datei erfolgreich generiert: %s\n', mat_filename);


% 1. Dummy-Ordner für 'functions' anlegen
func_dir = fullfile(scriptDir, 'functions');
if ~exist(func_dir, 'dir')
    mkdir(func_dir);
end

% 2. Temporäre Options-Datei für Modell 1 generieren
opt_file_m1 = 'temp_opts_m1.m';
fid = fopen(opt_file_m1, 'w');

fprintf(fid, 'function [modelname, paths, opts, prev_ident_pars] = temp_opts_m1()\n');
fprintf(fid, '    modelname = ''SF_Modell1.mat'';\n');

fprintf(fid, '    scriptDir = ''%s'';\n', scriptDir);
fprintf(fid, '    paths.models    = ''%s'';\n', model_path);
fprintf(fid, '    paths.results   = fullfile(scriptDir, ''results'');\n');
fprintf(fid, '    paths.functions = fullfile(scriptDir, ''functions'');\n');

fprintf(fid, '    opts.algorithm   = 1;\n');     %Algorithmus wählen -> Options schauen
fprintf(fid, '    opts.maxLietime  = 120;\n');   %Abbruchszeit
fprintf(fid, '    opts.decomp      = 0;\n');     %Zerlegung Hier nicht nötig weil kleines Modell
fprintf(fid, '    opts.checkObser  = 1;\n');     %Zusätzlich Beobachtbarkeit prüfen
fprintf(fid, '    opts.forcedecomp = 0;\n');     %Zerlegung auch wenn es das Modell nciht so sieht
fprintf(fid, '    opts.tol         = 1e-6;\n');  %Numerische Toleranz
fprintf(fid, '    opts.multiexp    = 0;\n');     %Gibt es Multiexperiment Daten -> Hier nur ein Batch also nein
fprintf(fid, '    opts.nnzDerU     = inf;\n');   %Wie oft dürfen die Stellgrüßen 0 sein
fprintf(fid, '    opts.nnzDerW     = 0;\n');     %Nöchmal für unbekannte Stellgrößen -> haben keine also 0
fprintf(fid, '    prev_ident_pars  = [];\n');    %Array für Parameter

fprintf(fid, 'end\n');
fclose(fid);

% 3. STRIKE-GOLDD Aufruf
try
    STRIKE_GOLDD(opt_file_m1);
catch ME
    fprintf('[FEHLER] Modell 1: %s\n', ME.message);
end

% 4. Aufräumen der Files
if exist(opt_file_m1, 'file')
    delete(opt_file_m1);
end


%% 2b. FISPO-Analyse: Modell 1 (Batch, Monod, mit O2)

fprintf('\n');
fprintf('=============================================================\n');
fprintf('  FISPO-Analyse: MODELL 1 — Batch, Monod, MIT O2\n');
fprintf('  Zustände  : cX, cGlc, cO2\n');
fprintf('  Parameter : mumax, KS, YXS, YXO, KLa, cO2stern\n');
fprintf('=============================================================\n\n');

% Modell aufstellen
syms cX cGlc cO2 real       
syms mumax KS YXS YXO KLa real  

x        = [cX; cGlc; cO2];                          % 1. Zustandsvektor (inkl. O2)
p        = [mumax; KS; YXS; YXO; KLa];     % 2. Parametervektor (inkl. O2-Parameter)
u        = [];                                       % 3. Bekannte Eingänge 
w        = [];                                       % 4. Unbekannte Eingänge / Störgrößen
ic       = [];                                       % 5. Anfangsbedingungen
known_ic = [];                                       % 6. Bekannte ICs

H = 1;
cO2stern = 100;

% Wachstumsrate und Differentialgleichungen (f)
mu = mumax * cGlc / (KS + cGlc);
f  = [ mu * cX;                                         % Biomasse
      -1/YXS * mu * cX;                                 % Glucose
       KLa * (cO2stern - cO2) - H * 1/YXO * mu * cX];       % O2-Bilanz

% Messgleichungen (h) - Annahme: Alle 3 Zustände können gemessen werden
h  = [cX; cGlc; cO2];              

% Als .mat speichern für Toolbox
mat_filename_m1_o2 = fullfile(model_path, 'SF_Modell1_O2.mat');
save(mat_filename_m1_o2, 'x', 'f', 'h', 'u', 'p', 'w', 'ic', 'known_ic');
fprintf('Modell-Datei erfolgreich generiert: %s\n', mat_filename_m1_o2);


% 1. Dummy-Ordner für 'functions' anlegen (zur Sicherheit)
func_dir = fullfile(scriptDir, 'functions');
if ~exist(func_dir, 'dir')
    mkdir(func_dir);
end

% 2. Temporäre Options-Datei für Modell 1 mit O2 generieren
opt_file_m1_o2 = 'temp_opts_m1_o2.m';
fid = fopen(opt_file_m1_o2, 'w');

fprintf(fid, 'function [modelname, paths, opts, prev_ident_pars] = temp_opts_m1_o2()\n');
fprintf(fid, '    modelname = ''SF_Modell1_O2.mat'';\n');

fprintf(fid, '    scriptDir = ''%s'';\n', scriptDir);
fprintf(fid, '    paths.models    = ''%s'';\n', model_path);
fprintf(fid, '    paths.results   = fullfile(scriptDir, ''results'');\n');
fprintf(fid, '    paths.functions = fullfile(scriptDir, ''functions'');\n');

fprintf(fid, '    opts.algorithm   = 1;\n');     
fprintf(fid, '    opts.maxLietime  = 120;\n');   
fprintf(fid, '    opts.decomp      = 0;\n');     
fprintf(fid, '    opts.checkObser  = 1;\n');     
fprintf(fid, '    opts.forcedecomp = 0;\n');     
fprintf(fid, '    opts.tol         = 1e-6;\n');  
fprintf(fid, '    opts.multiexp    = 0;\n');     
fprintf(fid, '    opts.nnzDerU     = inf;\n');   
fprintf(fid, '    opts.nnzDerW     = 0;\n');     
fprintf(fid, '    prev_ident_pars  = [];\n');    

fprintf(fid, 'end\n');
fclose(fid);

% 3. STRIKE-GOLDD Aufruf
try
    STRIKE_GOLDD(opt_file_m1_o2);
catch ME
    fprintf('[FEHLER] Modell 1 (mit O2): %s\n', ME.message);
end

% 4. Aufräumen der Files
if exist(opt_file_m1_o2, 'file')
    delete(opt_file_m1_o2);
end




%% 3. FISPO-Analyse: Modell 2 (Erweiterung Base/Ammonium, O2)

fprintf('\n');
fprintf('=============================================================\n');
fprintf('  FISPO-Analyse: MODELL 2 — Batch, Monod, O2 & pH-Regulation\n');
fprintf('  Zustände  : cX, cGlc, cAm, cBase, cO2\n');
fprintf('  Parameter : mumax, KS, YXS, YBam, YAmX, YXO, KLa, cO2stern\n');
fprintf('=============================================================\n\n');

% Modell aufstellen
syms cX cGlc cAm cBase cO2 real       
syms mumax KS YXS YBam YAmX YXO KLa cO2stern real  

x        = [cX; cGlc; cAm; cBase; cO2];                        % 1. Zustandsvektor (inkl. Ammonium, Base, O2)
p        = [mumax; KS; YXS; YBam; YAmX; YXO; KLa; cO2stern];   % 2. Parametervektor 
u        = [];                                                 % 3. Bekannte Eingänge 
w        = [];                                                 % 4. Unbekannte Eingänge / Störgrößen
ic       = [];                                                 % 5. Anfangsbedingungen
known_ic = [];                                                 % 6. Bekannte ICs

% Wachstumsrate und Differentialgleichungen (f) nach Modell 2 (Standard Monod)
mu = mumax * cGlc / (KS + cGlc);
f  = [ mu * cX;                                     % Biomasse
      - (1 / YXS) * mu * cX;                        % Glucose
      - YAmX * mu * cX;                             % Ammonium
        YBam * YAmX * mu * cX;                      % Base
        KLa * (cO2stern - cO2) - 1/YXO * mu * cX];  % O2

% Messgleichungen (h)
h  = [cX; cGlc; cAm; cBase; cO2];  



% Als .mat speichern für Toolbox
mat_filename_m2 = fullfile(model_path, 'SF_Modell2.mat');
save(mat_filename_m2, 'x', 'f', 'h', 'u', 'p', 'w', 'ic', 'known_ic');
fprintf('Modell-Datei erfolgreich generiert: %s\n', mat_filename_m2);


% 1. Dummy-Ordner für 'functions' anlegen (bereits in M1, aber zur Sicherheit)
func_dir = fullfile(scriptDir, 'functions');
if ~exist(func_dir, 'dir')
    mkdir(func_dir);
end

% 2. Temporäre Options-Datei für Modell 2 generieren
opt_file_m2 = 'temp_opts_m2.m';
fid = fopen(opt_file_m2, 'w');

fprintf(fid, 'function [modelname, paths, opts, prev_ident_pars] = temp_opts_m2()\n');
fprintf(fid, '    modelname = ''SF_Modell2.mat'';\n');

fprintf(fid, '    scriptDir = ''%s'';\n', scriptDir);
fprintf(fid, '    paths.models    = ''%s'';\n', model_path);
fprintf(fid, '    paths.results   = fullfile(scriptDir, ''results'');\n');
fprintf(fid, '    paths.functions = fullfile(scriptDir, ''functions'');\n');

fprintf(fid, '    opts.algorithm   = 1;\n');     
fprintf(fid, '    opts.maxLietime  = 120;\n');   
fprintf(fid, '    opts.decomp      = 0;\n');     
fprintf(fid, '    opts.checkObser  = 1;\n');     
fprintf(fid, '    opts.forcedecomp = 0;\n');     
fprintf(fid, '    opts.tol         = 1e-6;\n');  
fprintf(fid, '    opts.multiexp    = 0;\n');     
fprintf(fid, '    opts.nnzDerU     = inf;\n');   
fprintf(fid, '    opts.nnzDerW     = 0;\n');     
fprintf(fid, '    prev_ident_pars  = [];\n');    

fprintf(fid, 'end\n');
fclose(fid);

% 3. STRIKE-GOLDD Aufruf
try
    STRIKE_GOLDD(opt_file_m2);
catch ME
    fprintf('[FEHLER] Modell 2: %s\n', ME.message);
end

% 4. Aufräumen der Files
if exist(opt_file_m2, 'file')
    delete(opt_file_m2);
end




%% 4. FISPO-Analyse: Modell 3 ()





fprintf('\n=== Alle FISPO-Analysen abgeschlossen ===\n');
