%Executes the STRIKE-GOLDD toolbox for structural observability and identifiability (Task 3)
%
% main_03_fispo_Modell3_woEtOH
%
% FISPO-Analyse fuer Modell 3 — Untermodell OHNE Ethanol.
% Messgleichung h: Konzentrations-Form nach Modell3_mgl.m (ohne cEt).

clear; clc; close all;


%% 1. STRIKE-GOLDD Toolbox zum Pfad hinzufuegen

scriptDir  = pwd;
sg_path    = fullfile(fileparts(pwd), 'STRIKE-GOLDD');
model_path = fullfile(fileparts(pwd), 'Modelle');

if exist(sg_path, 'dir')
    addpath(genpath(sg_path));
    fprintf('[OK] STRIKE-GOLDD gefunden:\n     %s\n', sg_path);
else
    error(['[FEHLER] STRIKE-GOLDD nicht gefunden unter:\n  %s\n\n', ...
           'Erwartete Struktur: SPI_Projekt/STRIKE-GOLDD/strike_goldd.m\n'], sg_path);
end

if exist(model_path, 'dir')
    addpath(model_path);
    fprintf('[OK] Modell-Ordner gefunden:\n     %s\n', model_path);
else
    error('[FEHLER] Modell-Ordner nicht gefunden unter:\n  %s\n\n', model_path);
end


%% Bekannte (numerische) Konstanten
cO2_sat  = 7.14e-3;
DOTstern = 100;
cO2stern = DOTstern;
H        = 100 / cO2_sat;

cGlc_in = 450;
cAm_in  = 30;
cPh_in  = 24;


%% 2. FISPO-Analyse: Modell 3 — OHNE Ethanol

fprintf('\n');
fprintf('=============================================================\n');
fprintf('  FISPO-Analyse: MODELL 3 (woEtOH) — Fed-Batch ohne Ethanol\n');
fprintf('  Zustaende  : V, mX, mGlc, mAm, mPh, mB, DOT\n');
fprintf('  Parameter  : mumax, KS, YXS, YAmX, YPhX, YB_Am, KLa, YXO\n');
fprintf('  Messungen  : cX, cGlc, cNH4, cPO4, mB, V (mgl-Form)\n');
fprintf('=============================================================\n\n');

% --- Zustaende (kein mEt) ---
syms V mX mGlc mAm mPh mB DOT real

% --- Parameter ---
syms mumax KS YXS YAmX YPhX YB_Am KLa YXO real

% --- Bekannte Eingaenge (Feedraten) ---
syms uGlc uAm uPh uBase uAcid real

x = [V; mX; mGlc; mAm; mPh; mB; DOT];                   % 1. Zustandsvektor
p = [mumax; KS; YXS; YAmX; YPhX; YB_Am; KLa; YXO];      % 2. Parametervektor
u = [uGlc; uAm; uPh; uBase; uAcid];                     % 3. Bekannte Eingaenge
w = [];                                                 % 4. Unbekannte Eingaenge
ic       = [];                                          % 5. Anfangsbedingungen
known_ic = [];                                          % 6. Bekannte ICs

% --- Konzentrationen ---
cGlc = mGlc / V;
cAm  = mAm  / V;
cPh  = mPh  / V;

% --- Reaktionsrate (nur Wachstum auf Glucose, Gl. 16) ---
rX = mumax * cGlc / (KS + cGlc);

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
     V;        % Volumen
     DOT];      % DOT (wird gemessen -> macht KLa & YXO identifizierbar)

% Als .mat speichern
mat_filename = fullfile(model_path, 'SF_Modell3_woEtOH.mat');
save(mat_filename, 'x', 'f', 'h', 'u', 'p', 'w', 'ic', 'known_ic');
fprintf('Modell-Datei erfolgreich generiert: %s\n', mat_filename);


%% 3. Temporaere Options-Datei generieren und STRIKE-GOLDD aufrufen

opt_funname = 'temp_opts_SF_Modell3_woEtOH';
opt_file    = [opt_funname '.m'];

func_dir = fullfile(scriptDir, 'functions');
if ~exist(func_dir, 'dir')
    mkdir(func_dir);
end

fid = fopen(opt_file, 'w');
fprintf(fid, 'function [modelname, paths, opts, prev_ident_pars] = %s()\n', opt_funname);
fprintf(fid, '    modelname = ''SF_Modell3_woEtOH.mat'';\n');
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
fprintf(fid, '    opts.nnzDerU     = ones(1, %d);\n', numel(u)); % 1 Ableitung pro Feedrate (stueckweise konstant)
fprintf(fid, '    opts.nnzDerW     = [];\n');
fprintf(fid, '    prev_ident_pars  = [];\n');
fprintf(fid, 'end\n');
fclose(fid);

try
    STRIKE_GOLDD(opt_file);
catch ME
    fprintf('[FEHLER] Modell 3 (woEtOH): %s\n', ME.message);
end

if exist(opt_file, 'file')
    delete(opt_file);
end


fprintf('\n=== FISPO-Analyse (Modell 3 woEtOH) abgeschlossen ===\n');