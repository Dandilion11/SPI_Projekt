% OVP_Modell3.m
function ERGEBNIS = OVP_Modell3(VP, FM_old, invC)

% Indizes für endliche Schranken extrahieren
VP.CONS.xmin_idx = find(isfinite(VP.CONS.xmin));
VP.CONS.xmax_idx = find(isfinite(VP.CONS.xmax));

% Startprofil für die Optimierung (Glucose-Feedrate aus Zeile 6)
u0_search = VP.u(6, :);

% fmincon aufrufen
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 5000);
uOpt_glc = fmincon(@guete_ovp, u0_search, [], [], [], [], VP.CONS.umin, VP.CONS.umax, @cons_fcn, options);

% Optimales Profil im Ergebnis-Struct speichern
ERGEBNIS.t = VP.t;
ERGEBNIS.u = VP.u;
ERGEBNIS.u(6, :) = uOpt_glc;

% Abschließende Simulation mit optimalen Werten
[~, ERGEBNIS.x, ERGEBNIS.y] = simulation_fim_modell3(VP.t, VP.x0, ERGEBNIS.u, VP.p, invC);

    %% Verschachtelte Gütefunktion (A-Kriterium)
    function I = guete_ovp(u_current)
        % Aktuellen Vektor in die u-Matrix einsetzen
        VP.u(6, :) = u_current;
        
        % FIM für aktuellen Stellgrößenverlauf berechnen
        % HINWEIS: Funktion 'simulation_fim_modell3' muss im nächsten Schritt erstellt werden
        FM = simulation_fim_modell3(VP.t, VP.x0, VP.u, VP.p, invC) + FM_old;
        
        % Parameternormierung zur Vermeidung numerischer Schlechtkonditionierung
        FM_norm = diag(VP.p) * FM * diag(VP.p);
        
        % A-Kriterium: Minimierung der Spur der inversen FIM
        try
            invFM = inv(FM_norm);
            I = sum(diag(invFM));
        catch
            I = 1e6; % Strafterm bei singulärer Matrix
        end
    end

    %% Verschachtelte Nichtlineare Nebenbedingungen (Zustandsgrenzen)
    function [c, ceq] = cons_fcn(u_current)
        VP.u(6, :) = u_current;
        
        % Simulation ausführen, um Zustandstrajektorien zu erhalten
        [~, x_sim, ~] = simulation_fim_modell3(VP.t, VP.x0, VP.u, VP.p, invC);
        
        % Verletzung der minimalen und maximalen Zustandsschranken prüfen
        c = [ reshape(VP.CONS.xmin(VP.CONS.xmin_idx) * ones(1, length(VP.t)) - x_sim(VP.CONS.xmin_idx, :), [], 1); ...
              reshape(x_sim(VP.CONS.xmax_idx, :) - VP.CONS.xmax(VP.CONS.xmax_idx) * ones(1, length(VP.t)), [], 1) ];
        ceq = [];
    end

end