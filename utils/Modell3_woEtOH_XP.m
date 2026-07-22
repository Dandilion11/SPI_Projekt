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