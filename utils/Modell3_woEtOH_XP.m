function dxdt = Modell3_woEtOH_XP(t, x, u, p, DOTstern)
% Erweitertes System: Zustaende + Sensitivitaeten XP = dx/dtheta.
%   dx/dt  = f(x,u,p)
%   dXP/dt = df/dx * XP + df/dp        (Sensitivitaetsgleichungen)
%
% np = 10, weil das aktuelle Modell die zusaetzlichen Monod-Terme mit
% KAm und KPh enthaelt. Die im Fit fixierten Parameter bleiben hier drin
% -- ihre Spalten werden erst in der FIM herausgeschnitten.

nx = 7;  np = 10;

xs = x(1:nx);
XP = reshape(x(nx+1:end), nx, np);

dx = Modell3_woEtOH_10p(t, xs, u, p, DOTstern);
A  = Modell3_dfdx_10p(xs, p, DOTstern);      % nx x nx
B  = Modell3_dfdp_10p(xs, p, DOTstern);      % nx x np

dXP = A*XP + B;

dxdt = [dx(:); dXP(:)];
end

%{
function dxdt = Modell3_woEtOH_XP(t, x, u, p, DOTstern)
% Erweiterte DGL: Zustaende + Sensitivitaeten XP = dx/dp
nx = 7; np = 8;

x_state = x(1:nx);
XP      = reshape(x(nx+1:end), nx, np);

% Basis-DGL
dx_state = Modell3_woEtOH(t, x_state, u, p, DOTstern);

% Sensitivitaets-DGL: d/dt(XP) = dfdx*XP + dfdp
dfdx = Modell3_dfdx(x_state, p, DOTstern);
dfdp = Modell3_dfdp(x_state, p, DOTstern);
dXP  = dfdx * XP + dfdp;

dxdt = [dx_state; dXP(:)];
end
%}