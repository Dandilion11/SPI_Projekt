%% check_u_matrix_v2.m
%  Version 2. Aenderungen gegenueber v1:
%    - interp1-Absturz behoben: D.V wird sortiert und dedupliziert.
%      (Das ist selbst ein Befund: aufbereiten() sortiert Weight NICHT,
%       V0 = D.V.y(1) ist also evtl. gar nicht der Wert bei t_min.)
%    - NEU: u_Base gegen die gemessene Base-Zeitreihe m_B geprueft.
%      m_B ist definitionsgemaess das Integral von u_Base -- damit laesst
%      sich die Basenbilanz OHNE Fit testen.
%    - NEU: zeitvariable Zulaufkonzentrationen c_in werden gemeldet.
%    - NEU: Zeilenzahl von u wird geprueft (Def06/Def10 haben 11 statt 10).
%    - NEU: integrierte Feedmasse im Batch-Fenster von Modell 1/2
%      (Def03 t<=6 h, Def04 t<=8 h) -- dort wird "feedfrei" behauptet.
% =========================================================================
clear; clc;

kandidaten = { fullfile(pwd,'Daten','Daten_Processed','Processed_FedBatch_Modell3_MultiExp.mat'), ...
               fullfile(pwd,'..','Daten','Daten_Processed','Processed_FedBatch_Modell3_MultiExp.mat'), ...
               fullfile(pwd,'..','Matlab_Code','Daten','Daten_Processed','Processed_FedBatch_Modell3_MultiExp.mat') };
matfile = '';
for k = 1:numel(kandidaten)
    if exist(kandidaten{k},'file'), matfile = kandidaten{k}; break; end
end
if isempty(matfile), error('MAT-Datei nicht gefunden.'); end
fprintf('Datei: %s\n', matfile);
S = load(matfile);
Alle = [S.TrainSet, S.ValSet];

feedRows = [2 4 6 9 10];
concRows = [3 5 7];
concName = {'cAm_in','cPh_in','cGlc_in'};

%% ======================================================================
%  A. Zeilenzahl und zeitvariable Zulaufkonzentrationen
%  ======================================================================
fprintf('\n########## A. STRUKTUR VON u ##########\n');
for k = 1:numel(Alle)
    D = Alle(k);  u = D.u;
    fprintf('\n%-12s size(u) = %d x %d', D.name, size(u,1), size(u,2));
    if size(u,1) ~= 10
        fprintf('   [ACHTUNG] nicht 10 Zeilen!');
        for r = 11:size(u,1)
            fprintf('  (Zeile %d: nnz = %d, max = %.4g)', r, nnz(u(r,:)), max(u(r,:)));
        end
    end
    fprintf('\n');
    for i = 1:numel(concRows)
        r = concRows(i);
        vals = unique(u(r,:));
        if numel(vals) > 1
            fprintf('   [ACHTUNG] %s ist ZEITVARIABEL: %s\n', ...
                    concName{i}, mat2str(round(vals,3)));
            % Wann wechselt sie und was ist der Startwert?
            j = find(diff(u(r,:)) ~= 0);
            fprintf('              Startwert = %.4g | Wechsel bei t = %s h\n', ...
                    u(r,1), mat2str(round(u(1, j+1), 3)));
            fprintf(['              -> PRUEFEN: liest Modell3_woEtOH_10p c_in ', ...
                     'zeitabhaengig oder nur u(%d,1)?\n'], r);
        end
    end
end

%% ======================================================================
%  B. Weight-Zeitreihe: sortiert? dupliziert? stimmt V0?
%  ======================================================================
fprintf('\n\n########## B. WEIGHT / VOLUMEN ##########\n');
fprintf('%-12s %8s %8s %10s %10s %10s %8s %8s\n', ...
        'Exp','n(V)','sortiert','V.y(1)','V bei tmin','V(tmax)','dupl.','Reihe?');
for k = 1:numel(Alle)
    D = Alle(k);
    tV = D.V.t(:);  yV = D.V.y(:);
    sortiert = issorted(tV);
    nDup     = numel(tV) - numel(unique(tV));

    [tS, ord] = sort(tV);  yS = yV(ord);
    [tU, ia]  = unique(tS, 'stable');  yU = yS(ia);   % erster Wert je Zeit

    fprintf('%-12s %8d %8d %10.3f %10.3f %10.3f %8d %8d\n', ...
            D.name, numel(tV), sortiert, yV(1), yU(1), yU(end), nDup, D.V.isTimeSeries);
    if abs(yV(1) - yU(1)) > 1e-6
        fprintf(['   [FEHLER] V0 = D.V.y(1) = %.3f, aber der frueheste ', ...
                 'Zeitpunkt (t = %.3f h) hat V = %.3f.\n'], yV(1), tU(1), yU(1));
        fprintf('            x0 = [V0; c0*V0; ...] ist damit komplett falsch skaliert.\n');
    end
    % zur Weiterverwendung merken
    Alle(k).Vs.t = tU;  Alle(k).Vs.y = yU;
end

%% ======================================================================
%  C. Volumenbilanz (ZOH), jetzt mit bereinigter Weight-Reihe
%  ======================================================================
fprintf('\n\n########## C. VOLUMENBILANZ ##########\n');
fprintf('%-12s %8s %10s %10s %10s %10s %10s\n', ...
        'Exp','V0','Zufuhr','Proben','V_pred','V_mess','Abw. [%]');
for k = 1:numel(Alle)
    D  = Alle(k);  u = D.u;  tu = u(1,:);
    fr = feedRows(feedRows <= size(u,1));
    q  = sum(u(fr,:),1);
    dt = diff(tu);
    Vfeed = [0, cumsum(q(1:end-1).*dt)];
    tp = D.Probe.BatchAge(:);  vp = D.Probe.Volumen(:);
    Vprobe = arrayfun(@(tt) sum(vp(tp <= tt)), tu);
    V0 = D.x0(1);
    Vpred = V0 + Vfeed(end) - Vprobe(end);

    if D.V.isTimeSeries && numel(D.Vs.t) > 5
        tq = min(tu(end), max(D.Vs.t));
        Vmess = interp1(D.Vs.t, D.Vs.y, tq, 'linear', 'extrap');
        abw = 100*(Vpred - Vmess)/Vmess;
        fprintf('%-12s %8.3f %10.3f %10.3f %10.3f %10.3f %10.1f\n', ...
                D.name, V0, Vfeed(end), Vprobe(end), Vpred, Vmess, abw);
        if abs(abw) > 15
            fprintf('   [WARNUNG] >15 %% -> Zeilenbelegung, Einheiten oder V0 pruefen.\n');
        end
    else
        fprintf('%-12s %8.3f %10.3f %10.3f %10.3f %10s %10s\n', ...
                D.name, V0, Vfeed(end), Vprobe(end), Vpred, '--', '--');
    end
end

%% ======================================================================
%  D. Base: u_Base gegen die gemessene Zeitreihe m_B
%     m_B ist per Definition das Integral von u_Base. Stimmen beide
%     ueberein, ist m_B KEINE unabhaengige Messgroesse, sondern ein
%     bekannter Eingang -- und die Basenbilanz laesst sich ohne Fit pruefen.
%  ======================================================================
fprintf('\n\n########## D. BASE: EINGANG vs. MESSUNG ##########\n');
fprintf('%-12s %12s %12s %12s %10s\n', ...
        'Exp','int(u_Base)','m_B gemessen','m_B(0)','Verhaeltnis');
for k = 1:numel(Alle)
    D = Alle(k);  u = D.u;  tu = u(1,:);
    if size(u,1) < 10, fprintf('%-12s  (keine Zeile 10)\n', D.name); continue; end
    ub = u(10,:);
    dt = diff(tu);
    intU = sum(ub(1:end-1).*dt);

    mB = D.Base.y(:);
    dmB = mB(end) - mB(1);

    if intU < 1e-9 && dmB > 1e-3
        fprintf('%-12s %12.4f %12.4f %12.4f %10s   [WIDERSPRUCH] u_Base = 0, aber Base gemessen!\n', ...
                D.name, intU, dmB, mB(1), '--');
    else
        fprintf('%-12s %12.4f %12.4f %12.4f %10.3f\n', ...
                D.name, intU, dmB, mB(1), dmB/max(intU,eps));
    end
end
fprintf(['\nInterpretation: Verhaeltnis ~1 -> m_B = int(u_Base), d.h. die Base ist\n', ...
         'ein bekannter EINGANG. Dann laesst sich Y_B/Am*Y_Am/X direkt aus\n', ...
         'int(u_Base) / (erzeugte Biomasse) ablesen -- ohne Optimierung.\n']);

%% ======================================================================
%  E. Massenbilanz im Batch-Fenster (das eigentliche Def07-Raetsel)
%  ======================================================================
fprintf('\n\n########## E. MASSENBILANZ BATCH-FENSTER ##########\n');
fenster = struct('name',{'RamScDef03','RamScDef04','RamScDef06','RamScDef10','RamScDef07'}, ...
                 'tEnd',{6, 8, 8, 8, 8});
for f = 1:numel(fenster)
    k = find(strcmp({Alle.name}, fenster(f).name), 1);
    if isempty(k), continue; end
    D = Alle(k);  u = D.u;  tu = u(1,:);  tE = fenster(f).tEnd;

    fprintf('\n--- %s, Fenster t <= %.1f h ---\n', D.name, tE);
    kan = { 6, 7, D.Glucose,  'Glucose ';
            2, 3, D.Ammonium, 'Ammonium';
            4, 5, D.Phosphat, 'Phosphat' };
    for i = 1:size(kan,1)
        rU = kan{i,1};  rC = kan{i,2};  m = kan{i,3};  nm = kan{i,4};
        dt   = diff(tu);
        sel  = tu(1:end-1) < tE;
        dtc  = min(dt(sel), tE - tu([sel false]));
        mIn  = sum(u(rU,[sel false]) .* u(rC,[sel false]) .* dtc);

        ok = m.t <= tE;
        if nnz(ok) >= 2
            if D.V.isTimeSeries && numel(D.Vs.t) > 5
                Vm = interp1(D.Vs.t, D.Vs.y, m.t(ok), 'linear', 'extrap');
            else
                Vm = D.x0(1)*ones(nnz(ok),1);
            end
            mm = m.y(ok).*Vm(:);
            dMess = mm(end) - mm(1);
        else
            dMess = NaN;
        end
        fprintf('   %s: zugefuehrt = %8.2f g | Messung aendert sich um %+8.2f g', ...
                nm, mIn, dMess);
        if mIn > 1 && dMess < 0
            fprintf('   [WIDERSPRUCH]');
        elseif mIn > 1
            fprintf('   -> Verbrauch = %.2f g', mIn - dMess);
        end
        fprintf('\n');
    end
    % Biomasse zum Vergleich
    okX = D.Biomasse.t <= tE;
    if nnz(okX) >= 2
        if D.V.isTimeSeries && numel(D.Vs.t) > 5
            VmX = interp1(D.Vs.t, D.Vs.y, D.Biomasse.t(okX), 'linear', 'extrap');
        else
            VmX = D.x0(1)*ones(nnz(okX),1);
        end
        mX = D.Biomasse.y(okX).*VmX(:);
        fprintf('   Biomasse: +%.2f g gebildet -> impliziertes Y_XS = %.3f g/g\n', ...
                mX(end)-mX(1), NaN);
    end
end
fprintf(['\nLies Zeile "Verbrauch" gegen die gebildete Biomasse: Verbrauch/dX ist\n', ...
         '1/Y_XS. Kommt dabei etwas weit weg von 0.14 g/g heraus, stimmt entweder\n', ...
         'u oder die Messung in diesem Fenster nicht.\n']);

%% ======================================================================
%  F. Warnung fuer Modell 1 / 2
%  ======================================================================
fprintf('\n\n########## F. SIND DIE M1/M2-FENSTER FEEDFREI? ##########\n');
paare = {'RamScDef03', 6; 'RamScDef04', 8};
for f = 1:size(paare,1)
    k = find(strcmp({Alle.name}, paare{f,1}), 1);
    if isempty(k), continue; end
    D = Alle(k);  u = D.u;  tu = u(1,:);  tE = paare{f,2};
    aktiv = any(u([2 4 6],:) > 0, 1) & tu < tE;
    fprintf('%s (t <= %d h): %d von %d u-Spalten mit aktivem Feed', ...
            D.name, tE, nnz(aktiv), nnz(tu < tE));
    if any(aktiv)
        fprintf(' | erster Feed bei t = %.2f h  [Fenster ist NICHT feedfrei]\n', ...
                tu(find(aktiv,1)));
    else
        fprintf('  [feedfrei, ok]\n');
    end
end