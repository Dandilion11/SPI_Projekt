function X = sim_m3_sample(t,x0,u,p,DOTstern,Probe)

t = t(:);                 % Auswertezeiten als Spaltenvektor

tp = Probe.BatchAge(:);          % Probenzeitpunkte
vp = Probe.Volumen(:);          % jeweilige Probenvolumina
tf = u(1,:).';            % Zeitpunkte der Feed-Änderungen

t0 = u(1,1);              % Startzeit der Simulation

I = tp >= t0 & tp <= max(t);  % nur Proben im Simulationsbereich
tp = tp(I);
vp = vp(I);

tf = tf(tf >= t0 & tf <= max(t));  % nur Feed-Änderungen im Simulationsbereich

te = unique([t0; t; tf; tp]);      % alle Zeiten, an denen gestoppt werden muss

x = x0(:);                         % aktueller Zustand als Spaltenvektor
X = zeros(numel(t),numel(x0));      % Speicher für Simulationsergebnisse

opt = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

for k = 1:numel(te)

    if k > 1
        % Integration vom letzten Ereignis bis zum aktuellen Ereignis
        [~,Xs] = ode15s(@(tt,xx) Modell3_woEtOH(tt,xx,u,p,DOTstern), ...
                        te(k-1:k), x, opt);

        x = Xs(end,:).';           % Zustand am Ende des Intervalls übernehmen
    end

    % Zustand vor möglicher Probenentnahme speichern
    X(t == te(k),:) = repmat(x.',sum(t == te(k)),1);

    % falls an diesem Zeitpunkt Probe: Zustand sprunghaft ändern
    for j = find(tp == te(k)).'
        x = probe_m3(x,vp(j));
    end
end

end