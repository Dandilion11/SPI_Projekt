function X = sim_m3_sample_10p(t, x0, u, p, DOTstern, Probe)
% SIM_M3_SAMPLE_10P  Simuliert Modell 3 inklusive der Probenahme-Spruenge.
%
%   X = sim_m3_sample_10p(t, x0, u, p, DOTstern, Probe)
%   liefert die Zustaende zu den Zeitpunkten t (Zeilen = Zeitpunkte).
%
% Integriert wird segmentweise zwischen den echten Unstetigkeiten, also
% Feed-Sprung und Probenahme. Die AUSWERTEZEITEN sind bewusst KEINE
% Segmentgrenzen: sonst startet ode15s an jedem Messpunkt neu und das
% Ergebnis haengt vom Auswerteraster ab (Gitterabhaengigkeit 3.6e-1 statt
% 1e-8).

t  = t(:);
tp = Probe.BatchAge(:);   vp = Probe.Volumen(:);
tf = u(1,:).';            t0 = u(1,1);

% Nur Ereignisse im Simulationsbereich
I  = tp >= t0 & tp <= max(t);   tp = tp(I);   vp = vp(I);
tf = tf(tf >= t0 & tf <= max(t));
tb = unique([t0; tf; tp; max(t)]);        % Segmentgrenzen

x = x0(:);
X = zeros(numel(t), numel(x0));

% NonNegative haelt die Zustaende positiv (das Modell selbst klammert
% nicht). Das Schrittlimit bricht pathologische Parametersaetze ab, statt
% den Optimierer minutenlang haengen zu lassen -- deterministisch, im
% Gegensatz zu einem Zeitlimit.
opt = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'NonNegative',1:7, 'MaxStep',1.0, ...
             'OutputFcn', @(tt,yy,flag) ode_steplimit(tt,yy,flag,2e5));

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

    for i = find(sel).'
        X(i,:) = Xs(find(tout == t(i), 1), :);
    end

    x = Xs(end,:).';                      % Zustand am Segmentende
    for j = find(tp == te).'              % danach Probenahme
        x = probe_m3(x, vp(j));
    end
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