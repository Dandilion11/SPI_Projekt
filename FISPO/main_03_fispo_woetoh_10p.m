%Executes the STRIKE-GOLDD toolbox for structural observability and identifiability (Task 3)
%
% main_03_fispo_Modell3_woEtOH
%
% FISPO-Analyse fuer Modell 3 - Untermodell OHNE Ethanol.
% Messgleichung h: Konzentrations-Form nach Modell3_mgl.m (ohne cEt).
%
% Ordnerstruktur:
%   <Projekt>/Matlab_Code/FISPO/main_03_fispo_woetoh_10p.m   (dieses Skript)
%   <Projekt>/Matlab_Code/FISPO/results                      (Ergebnisse)
%   <Projekt>/Matlab_Code/STRIKE-GOLDD                       (Toolbox)
%   <Projekt>/Matlab_Code/Modelle                            (Modell-.mat)

clear; clc; close all;


%% 1. Pfade aufloesen (unabhaengig vom aktuellen MATLAB-Arbeitsordner)

scriptDir = fileparts(mfilename('fullpath'));   % .../Matlab_Code/FISPO
if isempty(scriptDir)                           % Fallback (z.B. Copy&Paste in Command Window)
    scriptDir = pwd;
end
code_root  = fileparts(scriptDir);              % .../Matlab_Code
sg_path    = fullfile(code_root, 'STRIKE-GOLDD');
model_path = fullfile(code_root, 'Modelle');

% WICHTIG: STRIKE-GOLDD legt die Zwischenergebnisse der Lie-Ableitungen in
% einem Ordner "results" RELATIV ZUM AKTUELLEN ARBEITSORDNER ab (nicht in
% paths.results). Deshalb wird hier in den Skript-Ordner gewechselt, sodass
% pwd/results == <FISPO>/results ist - der Ordner, der bereits existiert.
oldDir  = pwd;
cleanUp = onCleanup(@() cd(oldDir));            % Arbeitsordner am Ende zuruecksetzen
cd(scriptDir);

results_dir = fullfile(scriptDir, 'results');   % vorhandener Ergebnis-Ordner

if exist(sg_path, 'dir')
    addpath(genpath(sg_path));
    fprintf('[OK] STRIKE-GOLDD gefunden:\n     %s\n', sg_path);
else
    error(['[FEHLER] STRIKE-GOLDD nicht gefunden unter:\n  %s\n\n', ...
        'Erwartete Struktur: <projekt>/Matlab_Code/STRIKE-GOLDD/strike_goldd.m\n'], sg_path);
end

if exist(model_path, 'dir')
    addpath(model_path);
    fprintf('[OK] Modell-Ordner gefunden:\n     %s\n', model_path);
else
    error('[FEHLER] Modell-Ordner nicht gefunden unter:\n  %s\n', model_path);
end

if exist(results_dir, 'dir')
    fprintf('[OK] Ergebnis-Ordner gefunden:\n     %s\n', results_dir);
else
    error(['[FEHLER] Ergebnis-Ordner nicht gefunden unter:\n  %s\n\n', ...
        'STRIKE-GOLDD benoetigt diesen Ordner zum Speichern der Lie-Ableitungen.\n'], results_dir);
end

% paths.functions: Ordner mit den Toolbox-Hilfsfunktionen (kein neuer Ordner)
func_dir = fullfile(sg_path, 'functions');
if ~exist(func_dir, 'dir')
    func_dir = sg_path;
end


%% Bekannte (numerische) Konstanten
cO2_sat  = 7.14e-3;
DOTstern = 100;
cO2stern = DOTstern;
H        = 100 / cO2_sat;

cGlc_in = 450;
cAm_in  = 30;
cPh_in  = 24;


%% 2. FISPO-Analyse: Modell 3 - OHNE Ethanol

fprintf('\n');
fprintf('=============================================================\n');
fprintf('  FISPO-Analyse: MODELL 3 (woEtOH 10p) - Fed-Batch ohne Ethanol\n');
fprintf('  Zustaende  : V, mX, mGlc, mAm, mPh, mB, DOT\n');
fprintf('  Parameter  : mumax, KS, YXS, YAmX, YPhX, YB_Am, KLa, YXO, KAm, KPh\n');
fprintf('  Messungen  : cX, cGlc, cNH4, cPO4, mB, V, DOT (mgl-Form)\n');
fprintf('=============================================================\n\n');

% --- Zustaende (kein mEt) ---
syms V mX mGlc mAm mPh mB DOT real

% --- Parameter ---
syms mumax KS YXS YAmX YPhX YB_Am KLa YXO KAm KPh real

% --- Bekannte Eingaenge (Feedraten) ---
syms uGlc uAm uPh uBase uAcid real

x = [V; mX; mGlc; mAm; mPh; mB; DOT];                        % 1. Zustandsvektor
p = [mumax; KS; YXS; YAmX; YPhX; YB_Am; KLa; YXO; KAm; KPh]; % 2. Parametervektor
u = [uGlc; uAm; uPh; uBase; uAcid];                          % 3. Bekannte Eingaenge
w = [];                                                      % 4. Unbekannte Eingaenge
ic       = [];                                               % 5. Anfangsbedingungen
known_ic = [];                                               % 6. Bekannte ICs

% --- Konzentrationen ---
cGlc = mGlc / V;
cAm  = mAm  / V;
cPh  = mPh  / V;

% --- Reaktionsrate (nur Wachstum auf Glucose, Gl. 16) ---
rX = mumax * (cGlc/(KS + cGlc)) * (cAm/(KAm + cAm)) * (cPh/(KPh + cPh));

% --- Massenbilanzen (uout = 0) ---
f = [ uGlc + uAm + uPh + uBase + uAcid;                 % V
      rX * mX;                                          % mX
      (-1/YXS*rX)*mX + cGlc_in*uGlc;                    % mGlc
      -YAmX*rX*mX + cAm_in*uAm;                         % mAm
      -YPhX*rX*mX + cPh_in*uPh;                         % mPh
      YB_Am*YAmX*rX*mX;                                 % mB
      KLa*(cO2stern - DOT) - H*(1/YXO)*rX*(mX/V) ];     % DOT

% --- Messgleichungen (h): Konzentrations-Form nach mgl, ohne cEt ---
h = [mX/V;      % cX
     mGlc/V;    % cGlc
     mAm/V;     % cNH4
     mPh/V;     % cPO4
     mB;        % mB
     V;         % Volumen
     DOT];      % DOT (wird gemessen -> macht KLa & YXO identifizierbar)

% Als .mat speichern
mat_filename = fullfile(model_path, 'SF_Modell3_woEtOH_10p.mat');
save(mat_filename, 'x', 'f', 'h', 'u', 'p', 'w', 'ic', 'known_ic');
fprintf('Modell-Datei erfolgreich generiert: %s\n', mat_filename);


%% 3. Temporaere Options-Datei generieren und STRIKE-GOLDD aufrufen

opt_funname = 'temp_opts_SF_Modell3_woEtOH_10p';
opt_file    = [opt_funname '.m'];               % wird in scriptDir angelegt

% Backslashes maskieren, damit fprintf sie nicht als Escape-Sequenz liest
esc = @(s) strrep(s, '\', '\\');

fid = fopen(opt_file, 'w');
fprintf(fid, 'function [modelname, paths, opts, prev_ident_pars] = %s()\n', opt_funname);
fprintf(fid, '    modelname = ''SF_Modell3_woEtOH_10p.mat'';\n');
fprintf(fid, '    paths.models    = ''%s'';\n', esc(model_path));
fprintf(fid, '    paths.results   = ''%s'';\n', esc(results_dir));
fprintf(fid, '    paths.functions = ''%s'';\n', esc(func_dir));
fprintf(fid, '    opts.algorithm   = 1;\n');
fprintf(fid, '    opts.maxLietime  = 600;\n');   % war 120
fprintf(fid, '    opts.decomp      = 0;\n');
fprintf(fid, '    opts.checkObser  = 1;\n');
fprintf(fid, '    opts.forcedecomp = 0;\n');
fprintf(fid, '    opts.tol         = 1e-6;\n');
fprintf(fid, '    opts.multiexp    = 0;\n');
fprintf(fid, '    opts.nnzDerU     = ones(1, %d);\n', numel(u)); % 1 Ableitung pro Feedrate (stueckweise konstant)
fprintf(fid, '    opts.nnzDerW     = [];\n');
fprintf(fid, '    prev_ident_pars  = [];\n');
fprintf(fid, 'end\n');
fclose(fid);

addpath(scriptDir);
rehash;   % neue Options-Funktion sofort sichtbar machen

try
    STRIKE_GOLDD(opt_file);
catch ME
    fprintf('[FEHLER] Modell 3 (woEtOH 10p): %s\n', ME.message);
end

if exist(opt_file, 'file')
    delete(opt_file);
end

fprintf('\nErgebnisse liegen in:\n     %s\n', results_dir);
fprintf('\n=== FISPO-Analyse (Modell 3 woEtOH 10p) abgeschlossen ===\n');
