function X = sim_m3_sample_10p(t, x0, u, p, DOTstern, Probe, profil)
% SIM_M3_SAMPLE_10P  Simuliert Modell 3 inklusive der Probenahme-Spruenge.
%
%   X = sim_m3_sample_10p(t, x0, u, p, DOTstern, Probe)
%   X = sim_m3_sample_10p(..., profil)      profil = 'exact'|'fast'|'coarse'
%
% Integriert wird segmentweise zwischen den Unstetigkeiten, also Feed-Sprung
% und Probenahme. Die AUSWERTEZEITEN sind bewusst KEINE Segmentgrenzen:
% sonst startet ode15s an jedem Messpunkt neu und das Ergebnis haengt vom
% Auswerteraster ab (Gitterabhaengigkeit 3.6e-1 statt 1e-8).
%
% BESCHLEUNIGUNGEN gegenueber der Vorfassung:
%
%   1. SEGMENTAUSDUENNUNG (groesster Hebel). Die Feeds werden von
%      PI-geregelten Pumpen gestellt und als RAMPE mitgeschrieben -- fast
%      jede Spalte von u unterscheidet sich minimal von der vorigen. Bisher
%      wurde daraus je ein ode15s-Neustart (284 bzw. 361 pro Lauf), obwohl
%      es keine echten Spruenge sind. Jetzt wird nur noch dort neu gestartet,
%      wo sich eine Stellgroesse um mehr als du_rel ihres Maximums aendert.
%      WICHTIG: das MODELL bekommt weiterhin die vollstaendige u-Matrix, die
%      zugefuehrten Mengen aendern sich also nicht. Ausgeduennt wird nur die
%      Segmentierung -- eine kontrollierte Genauigkeitseinbusse, die gegen
%      profil='exact' zu pruefen ist.
%
%   2. u wird verlustfrei auf Spalten mit tatsaechlicher Aenderung reduziert
%      (exakte Duplikate). Verkuerzt zusaetzlich die Suche
%      find(t >= u(1,:),1,'last') in jedem RHS-Aufruf.
%
%   3. Analytische Jacobimatrix fuer ode15s aus Modell3_dfdx_10p. Ohne sie
%      baut der Solver df/dx bei jeder Neufaktorisierung numerisch auf
%      (7 zusaetzliche RHS-Auswertungen). Signatur wird beim ersten Aufruf
%      ermittelt und geprueft; schlaegt das fehl, laeuft alles wie bisher.
%
%   4. Toleranzprofile, siehe build_opts.
%
% PROFIL SETZEN -- zwei Wege:
%
%   (a) Empfohlen, in main_04 die SIMFUN-Zeile aendern:
%           SIMFUN = @(tt,xx,uu,pp,dd,PP) ...
%                    sim_m3_sample_10p(tt,xx,uu,pp,dd,PP,'fast');
%       Funktioniert auch mit parfor / UseParallel.
%
%   (b) Ohne jede Aenderung am Aufrufer, ueber eine globale Variable:
%           global SIM_M3_PROFIL;  SIM_M3_PROFIL = 'coarse';
%       Wird NICHT an parallele Worker uebertragen.
%
%   Ohne Angabe gilt 'fast'. Das 7. Argument hat Vorrang vor der Globalen.

global SIM_M3_PROFIL %#ok<GVMIS>
persistent OPTCACHE JACFUN JACSTATE
if isempty(JACSTATE), JACSTATE = 0; end     % 0 ungeprueft, 1 ok, -1 aus

%% Profil bestimmen -------------------------------------------------------
if nargin < 7 || isempty(profil)
    if isempty(SIM_M3_PROFIL), profil = 'fast'; else, profil = SIM_M3_PROFIL; end
end
feld = lower(profil);
if isempty(OPTCACHE) || ~isfield(OPTCACHE, feld)
    OPTCACHE.(feld) = build_opts(feld);
end
cfg = OPTCACHE.(feld);

%% u verlustfrei verdichten (Modell sieht weiterhin jeden Feedwert) -------
u = compress_u(u);

%% Jacobimatrix (einmalig pruefen) ----------------------------------------
if JACSTATE == 0
    [JACFUN, JACSTATE] = probe_jacobian(x0, u, p, DOTstern);
end

opt = cfg.ode;                               % Struktur direkt erweitern,
if JACSTATE == 1                             % odeset je Aufruf waere teuer
    opt.Jacobian = @(tt,xx) JACFUN(tt, xx, u, p, DOTstern);
end

%% Segmentgrenzen ---------------------------------------------------------
t  = t(:);
tp = Probe.BatchAge(:);   vp = Probe.Volumen(:);
t0 = u(1,1);

% Nur SIGNIFIKANTE Feedwechsel werden Segmentgrenze
tf = signif_switch_times(u, cfg.du_rel);

I  = tp >= t0 & tp <= max(t);   tp = tp(I);   vp = vp(I);
tf = tf(tf >= t0 & tf <= max(t));
tb = unique([t0; tf; tp; max(t)]);

x = x0(:);
X = zeros(numel(t), numel(x0));

% Startzeitpunkt: Zustand VOR einer eventuellen Probe bei t0
X(t == t0,:) = repmat(x.', nnz(t == t0), 1);
for j = find(tp == t0).'
    x = probe_m3(x, vp(j));
end

for k = 2:numel(tb)
    ta = tb(k-1);   te = tb(k);
    sel  = t > ta & t <= te;              % Auswertezeiten in (ta, te]
    tout = unique([ta; t(sel); te]);

    [~, Xs] = ode15s(@(tt,xx) Modell3_woEtOH_10p(tt,xx,u,p,DOTstern), ...
                     tout, x, opt);

    % Bei genau zwei Zeitpunkten gibt ode15s alle internen Schritte zurueck
    if numel(tout) == 2
        Xs = Xs([1 end],:);
    end
    if size(Xs,1) ~= numel(tout)
        return                            % abgebrochen -> X bleibt zu kurz
    end

    idxT = find(sel);                     % Zuordnung in einem Rutsch
    if ~isempty(idxT)
        [~, loc] = ismember(t(idxT), tout);
        X(idxT,:) = Xs(loc,:);
    end

    x = Xs(end,:).';                      % Zustand am Segmentende
    for j = find(tp == te).'               % danach Probenahme
        x = probe_m3(x, vp(j));
    end
end
end


%% ======================================================================
function cfg = build_opts(profil)
% Toleranzprofile. RelTol/AbsTol und der Differenzenschritt der PI haengen
% zusammen: ein kleinerer FD-Schritt verlangt eine genauere Integration.
% du_rel steuert die Segmentausduennung (0 = jede Aenderung ist Grenze).
switch profil
    case 'exact'      % identisch zur Vorfassung, fuer den Vergleichslauf
        rt = 1e-8;  at = 1e-10;  ms = 1.0;  nmax = 2e5;  du = 0;
    case 'fast'       % Standard fuer fmincon
        rt = 1e-6;  at = 1e-8;   ms = 5.0;  nmax = 1e5;  du = 0.02;
    case 'coarse'     % LHS-Screening, Multistart-Vorlauf
        rt = 1e-4;  at = 1e-6;   ms = 5.0;  nmax = 5e4;  du = 0.05;
    otherwise
        error('Unbekanntes Profil "%s" (exact|fast|coarse).', profil);
end

% NonNegative haelt die Zustaende positiv (das Modell klammert nicht). Das
% Schrittlimit bricht pathologische Parametersaetze ab, statt den Optimierer
% haengen zu lassen -- deterministisch, im Gegensatz zu einem Zeitlimit.
cfg.ode = odeset('RelTol',rt, 'AbsTol',at, 'NonNegative',1:7, 'MaxStep',ms, ...
                 'OutputFcn', @(tt,yy,flag) ode_steplimit(tt,yy,flag,nmax));
cfg.du_rel = du;
end


%% ======================================================================
function uc = compress_u(u)
% Verlustfrei: entfernt nur Spalten, die mit der vorigen identisch sind.
% Der Nullter-Ordnung-Halt bleibt exakt erhalten.
if size(u,2) < 2
    uc = u;  return;
end
d  = any(diff(u(2:end,:), 1, 2) ~= 0, 1);
uc = u(:, [true, d]);
end


%% ======================================================================
function tf = signif_switch_times(u, du_rel)
% Zeitpunkte, an denen sich eine Stellgroesse um mehr als du_rel ihres
% Maximums aendert. Bei du_rel = 0 sind das alle Spaltenzeitpunkte.
tf = u(1,:).';
if du_rel <= 0 || size(u,2) < 2
    return
end
U = u(2:end,:);
skal = max(abs(U), [], 2);            % Bezugsgroesse je Stellgroesse
skal(skal == 0) = 1;
rel  = max(abs(diff(U,1,2)) ./ skal, [], 1);
tf   = [u(1,1); u(1, [false, rel > du_rel]).'];
end


%% ======================================================================
function [fh, state] = probe_jacobian(x0, u, p, DOTstern)
% Ermittelt die Aufrufsignatur von Modell3_dfdx_10p und prueft das Ergebnis
% gegen eine numerische Jacobimatrix. Nur beim ersten Aufruf.
fh = [];  state = -1;
if exist('Modell3_dfdx_10p','file') ~= 2
    return
end

kandidaten = { @(tt,xx,uu,pp,dd) Modell3_dfdx_10p(tt,xx,uu,pp,dd), ...
               @(tt,xx,uu,pp,dd) Modell3_dfdx_10p(tt,xx,uu,pp),    ...
               @(tt,xx,uu,pp,dd) Modell3_dfdx_10p(xx,uu,pp,dd),    ...
               @(tt,xx,uu,pp,dd) Modell3_dfdx_10p(xx,pp) };

t0 = u(1,1);
x  = x0(:);
Jn = num_jacobian(t0, x, u, p, DOTstern);
skal = max(norm(Jn,1), 1);

for i = 1:numel(kandidaten)
    try
        Ja = kandidaten{i}(t0, x, u, p, DOTstern);
    catch
        continue
    end
    if ~isequal(size(Ja), [7 7]) || any(~isfinite(Ja(:)))
        continue
    end
    if norm(Ja - Jn, 1) / skal < 1e-4       % lockere Schranke: Jn ist selbst
        fh = kandidaten{i};                 % nur eine Naeherung
        state = 1;
        return
    end
end

warning(['sim_m3_sample_10p: Modell3_dfdx_10p konnte nicht als ', ...
         'ode15s-Jacobimatrix verwendet werden (Signatur oder Wert passt ', ...
         'nicht). Es wird ohne analytische Jacobimatrix gerechnet.']);
end


%% ======================================================================
function J = num_jacobian(t, x, u, p, DOTstern)
% Referenz-Jacobimatrix per zentraler Differenz, nur fuer die Pruefung.
n = numel(x);
J = zeros(n);
for j = 1:n
    h  = 1e-6 * max(abs(x(j)), 1);
    xp = x;  xp(j) = xp(j) + h;
    xm = x;  xm(j) = xm(j) - h;
    J(:,j) = (Modell3_woEtOH_10p(t,xp,u,p,DOTstern) - ...
              Modell3_woEtOH_10p(t,xm,u,p,DOTstern)) / (2*h);
end
end


%% ======================================================================
function status = ode_steplimit(tt, ~, flag, nmax)
% Bricht ab, wenn der Solver mehr als nmax Schritte braucht. Die
% Ausgabematrix ist dann zu kurz und der Aufrufer setzt J = 1e8.
persistent n
switch flag
    case 'init', n = 0;                 status = 0;
    case 'done',                        status = 0;
    otherwise,   n = n + numel(tt);     status = double(n > nmax);
end
end