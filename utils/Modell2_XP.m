function dxdt_ext = Modell2_XP(t, x_ext, p, kinetic)
% Kopplung von System-DGL und Sensitivitäts-DGL für Modell 2
% x_ext enthält [Zustände (5x1); Sensitivitäten XP (40x1)]

nx = 5; 
np = 8;

% Systemzustände und Sensitivitätsmatrix trennen
x  = x_ext(1:nx);
XP = reshape(x_ext(nx+1:end), nx, np);

% Zustandsschranke
x = max(x, 0);

% Auswertung der Einzelkomponenten
dxdt = Modell2(t, x, p, kinetic);
dfdx = Modell2_dfdx(x, p);
dfdp = Modell2_dfdp(x, p);

% Sensitivitäts-DGL: dXP/dt = df/dx * XP + df/dp
dXPdt = dfdx * XP + dfdp;

% Vektorisierte Rückgabe für ODE-Solver
dxdt_ext = [dxdt; dXPdt(:)];
end