function dxdt_ext = Modell2_XP_Monod(t, x_ext, p, DOTstern)
% Erweiterte DGL: Zustände + Sensitivitätsmatrix XP
% x_ext = [cX; cGlc; XP(:)] mit XP


nx = 5;
np = numel(p);
% Zustaende und Sensitivitäten trennen
x  = x_ext(1:nx);
XP = reshape(x_ext(nx+1:end), nx, np);

% Grenzen
x = max(x, 0);

% Modellgleichung
dxdt = Modell2(t, x, p,3, DOTstern);

% Jacobi-Matrizen
dfdx = Modell2_dfdx(x, p); 
dfdp = Modell2_dfdp(x, p, DOTstern);


% Sensitivitäts-DGL: dXP/dt = df/dx * XP + df/dp wie in Skript S. 34 Gl.
% 2.43
dXPdt = dfdx * XP + dfdp;

% Ergebnis zusammenbauen
dxdt_ext = [dxdt; dXPdt(:)];
end