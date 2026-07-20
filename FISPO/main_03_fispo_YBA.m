%Executes the STRIKE-GOLDD toolbox for structural observability and identifiability (Task 3)
%
% main_03_fispo_Modell3_YBA
%
% FISPO-Analyse fuer Modell 3 — Untermodell MIT Ethanol, aber mit
% gelumptem Base-Parameter YBA in der mB-Bilanz (statt YB_Am*YAmX).
% Messgleichung h: Konzentrations-Form nach Modell3_mgl.m.
%
% HINWEIS: In Modell3_YBA.m wird in der mAm-Bilanz YAmX verwendet, das aber
% NICHT im dortigen Parametervektor p steht. Damit das Modell vollstaendig
% und lauffaehig ist, wurde YAmX hier in p aufgenommen (p hat 14 Eintraege).
% Falls YAmX bewusst entfallen soll, bitte melden.

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


%% 2. FISPO-Analyse: Modell 3 — YBA (gelumpter Base-Parameter)

fprintf('\n');
fprintf('=============================================================\n');
fprintf('  FISPO-Analyse: MODELL 3 (YBA) — Fed-Batch mit Ethanol\n');
fprintf('  mB-Bilanz mit gelumptem Parameter YBA (= YB_Am*YAmX)\n');
fprintf('  Zustaende  : V, mX, mGlc, mAm, mPh, mB, DOT, mEt\n');
fprintf('  Parameter  : mumax, KS, YXS, YAmX, YPhX, YBA,\n');
fprintf('               mumax_EtP, mumax_EtX, YGlc_Et, YEt_X,\n');
fprintf('               KEt, KGlc_Et, KLa, YXO\n');
fprintf('  Messungen  : cX, cGlc, cNH4, cPO4, cEt, mB, V (mgl-Form)\n');
fprintf('=============================================================\n\n');

% --- Zustaende ---
syms V mX mGlc mAm mPh mB DOT mEt real

% --- Parameter (YBA statt YB_Am; YAmX ergaenzt, s. Hinweis oben) ---
syms mumax KS YXS YAmX YPhX YBA real
syms mumax_EtP mumax_EtX YGlc_Et YEt_X KEt KGlc_Et KLa YXO real

% --- Bekannte Eingaenge (Feedraten) ---
syms uGlc uAm uPh uBase uAcid real

x = [V; mX; mGlc; mAm; mPh; mB; DOT; mEt];              % 1. Zustandsvektor
p = [mumax; KS; YXS; YAmX; YPhX; YBA; ...             % 2. Parametervektor
     mumax_EtP; mumax_EtX; YGlc_Et; YEt_X; KEt; KGlc_Et; KLa; YXO];
u = [uGlc; uAm; uPh; uBase; uAcid];                     % 3. Bekannte Eingaenge
w = [];                                                 % 4. Unbekannte Eingaenge
ic       = [];                                          % 5. Anfangsbedingungen
known_ic = [];                                          % 6. Bekannte ICs

% --- Konzentrationen ---
cGlc = mGlc / V;
cEt  = mEt  / V;
cAm  = mAm  / V;
cPh  = mPh  / V;

% --- Reaktionsraten (Kraemer 2017, Gl. 16-18) ---
rX    = mumax     * cGlc / (KS + cGlc);
rEt_P = mumax_EtP * cGlc / (cGlc + KS);
rEt_X = mumax_EtX * (cEt / (cEt + KEt)) * (KGlc_Et / (cGlc + KGlc_Et));

% --- Massenbilanzen (uout = 0); mB-Bilanz mit gelumptem YBA ---
f = [ uGlc + uAm + uPh + uBase + uAcid;                 % V
      (rX + rEt_X) * mX;                                % mX
      (-1/YXS*rX - YGlc_Et*rEt_P)*mX + cGlc_in*uGlc;    % mGlc
      -YAmX*(rX+rEt_X)*mX + cAm_in*uAm;                 % mAm
      -YPhX*(rX+rEt_X)*mX + cPh_in*uPh;                 % mPh
      YBA*(rX+rEt_X)*mX;                                % mB  (YB_Am*YAmX -> YBA)
      KLa*(cO2stern - DOT) - H*(1/YXO)*(rX+rEt_X)*(mX/V); % DOT
      (rEt_P - YEt_X*rEt_X)*mX ];                       % mEt

% --- Messgleichungen (h): Konzentrations-Form nach mgl ---
h = [mX/V;      % cX
     mGlc/V;    % cGlc
     mAm/V;     % cNH4
     mPh/V;     % cPO4
     mEt/V;     % cEt
     mB;        % mB
     V;        % Volumen
     DOT];      % DOT (wird gemessen -> macht KLa & YXO identifizierbar)

% Als .mat speichern
mat_filename = fullfile(model_path, 'SF_Modell3_YBA.mat');
save(mat_filename, 'x', 'f', 'h', 'u', 'p', 'w', 'ic', 'known_ic');
fprintf('Modell-Datei erfolgreich generiert: %s\n', mat_filename);


%% 3. Temporaere Options-Datei generieren und STRIKE-GOLDD aufrufen

opt_funname = 'temp_opts_SF_Modell3_YBA';
opt_file    = [opt_funname '.m'];

func_dir = fullfile(scriptDir, 'functions');
if ~exist(func_dir, 'dir')
    mkdir(func_dir);
end

fid = fopen(opt_file, 'w');
fprintf(fid, 'function [modelname, paths, opts, prev_ident_pars] = %s()\n', opt_funname);
fprintf(fid, '    modelname = ''SF_Modell3_YBA.mat'';\n');
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
    fprintf('[FEHLER] Modell 3 (YBA): %s\n', ME.message);
end

if exist(opt_file, 'file')
    delete(opt_file);
end


fprintf('\n=== FISPO-Analyse (Modell 3 YBA) abgeschlossen ===\n');