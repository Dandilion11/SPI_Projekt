function dxdt_ext = Modell1_XP(t, x_ext, p, kinetic, withOxygen)
% Erweiterte DGL: Zustände + Sensitivitätsmatrix XP
% x_ext = [cX; cGlc; XP(:)] mit XP

nx = 2; np = 3;

% Zustaende und Sensitivitäten trennen
x  = x_ext(1:nx);
XP = reshape(x_ext(nx+1:end), nx, np);

% Grenzen
x = max(x, 0);

% Modellgleichung
dxdt = Modell1(t, x, p, kinetic, withOxygen);

% Jacobi-Matrizen
dfdx = Modell1_dfdx(x, p);
dfdp = Modell1_dfdp(x, p);

% Sensitivitäts-DGL: dXP/dt = df/dx * XP + df/dp wie in Skript S. 34 Gl.
% 2.43
dXPdt = dfdx * XP + dfdp;

% Ergebnis zusammenbauen
dxdt_ext = [dxdt; dXPdt(:)];
end