% OVP_Modell3.m
function ERGEBNIS = OVP_Modell3(VP, FM_old, invC)

VP.CONS.xmin_idx = find(isfinite(VP.CONS.xmin));
VP.CONS.xmax_idx = find(isfinite(VP.CONS.xmax));

u0_search = VP.u(6, :);   % Glucose-Feedrate als Optimierungsvariable

options = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
                       'MaxFunctionEvaluations',5000);
uOpt_glc = fmincon(@guete_ovp, u0_search, [], [], [], [], ...
                   VP.CONS.umin, VP.CONS.umax, @cons_fcn, options);

ERGEBNIS.t = VP.t;
ERGEBNIS.u = VP.u;
ERGEBNIS.u(6,:) = uOpt_glc;
[~, ERGEBNIS.x, ERGEBNIS.y] = ...
simulation_fim_modell3(VP.t, VP.x0, ERGEBNIS.u, VP.p, invC, VP.DOTstern);

    %% A-Kriterium
    function I = guete_ovp(u_current)
        VP.u(6,:) = u_current;
        FM = simulation_fim_modell3(VP.t, VP.x0, VP.u, VP.p, invC, VP.DOTstern) + FM_old;
        FM_norm = diag(VP.p) * FM * diag(VP.p);   % Normierung
        try
            I = sum(diag(inv(FM_norm)));
        catch
            I = 1e6;
        end
    end

    %% Zustandsschranken
    function [c, ceq] = cons_fcn(u_current)
        VP.u(6,:) = u_current;
        [~, x_sim, ~] = simulation_fim_modell3(VP.t, VP.x0, VP.u, VP.p, invC, VP.DOTstern);
        c = [ reshape(VP.CONS.xmin(VP.CONS.xmin_idx)*ones(1,length(VP.t)) - x_sim(VP.CONS.xmin_idx,:), [], 1); ...
              reshape(x_sim(VP.CONS.xmax_idx,:) - VP.CONS.xmax(VP.CONS.xmax_idx)*ones(1,length(VP.t)), [], 1) ];
        ceq = [];
    end
end
