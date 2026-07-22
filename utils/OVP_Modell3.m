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
%% FIM von Modell 3
function [FM, x_out, y_out] = simulation_fim_modell3(t, x0, u, p, invC, DOTstern)
% Simulation mit Sensitivitaeten und Berechnung der Fisher-Informationsmatrix
% invC : 6x6 (inverse Messkovarianz, 6 Messgroessen)
nx = 7; np = 8; ny = 6;

if nargin < 6 || isempty(DOTstern), DOTstern = 100; end

% Anfangsbedingungen (Zustaende + Nullsensitivitaeten)
XP0    = zeros(nx, np);          % 7 x 8
x0_ext = [x0(:); XP0(:)];        % 7 + 56 = 63 Elemente

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~, X_ext] = ode15s(@(tt,xx) Modell3_woEtOH_XP(tt,xx,u,p,DOTstern), t, x0_ext, opts);

FM    = zeros(np, np);
x_out = zeros(nx, numel(t));
y_out = zeros(ny, numel(t));

for k = 1:numel(t)
    xk         = X_ext(k, 1:nx).';
    x_out(:,k) = xk;
    y_out(:,k) = Modell3_mgl(xk);

    XP_k     = reshape(X_ext(k, nx+1:end), nx, np);
    dhdx     = Modell3_dmgldx(xk);
    dhdtheta = dhdx * XP_k;              % dhdp = 0
    FM       = FM + dhdtheta' * invC * dhdtheta;
end
end
