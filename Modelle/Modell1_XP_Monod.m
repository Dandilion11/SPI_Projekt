function dxdt_ext = Modell1_XP_Monod(t,nx, np, x_ext, p, withOxygen, KLaConst)
% Erweiterte DGL: Zustände + Sensitivitätsmatrix XP
% x_ext = [cX; cGlc; XP(:)] mit XP


% Zustaende und Sensitivitäten trennen
x  = x_ext(1:nx);
XP = reshape(x_ext(nx+1:end), nx, np);

% Grenzen
x = max(x, 0);

% Modellgleichung
dxdt = Modell1(t, x, p,3, withOxygen, KLaConst);

% Jacobi-Matrizen
if withOxygen
    dfdx = Modell1_dfdx_MonodO2(x, p); 
    dfdp = Modell1_dfdp_MonodO2(x, p);
else 
    dfdx = Modell1_dfdx_Monod(x, p);
    dfdp = Modell1_dfdp_Monod(x, p);
end


% Sensitivitäts-DGL: dXP/dt = df/dx * XP + df/dp wie in Skript S. 34 Gl.
% 2.43
dXPdt = dfdx * XP + dfdp;

% Ergebnis zusammenbauen
dxdt_ext = [dxdt; dXPdt(:)];
end