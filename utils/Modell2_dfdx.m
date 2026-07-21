function dfdx = Modell2_dfdx(x, p)
% Jacobimatrix df/dx für Modell2 (Monod-Kinetik)
% Zustände: x = [cX, cGlc, cAm, cBase, cO2]
% Parameter: p = [mumax, KS, YXS, YBam, YAmX, YXO, KLa]

cO2_sat = 7.14e-3;
H       = 100 / cO2_sat;

cX    = x(1);
cGlc  = x(2);
cO2   = x(5);

mumax = p(1);
KS    = p(2);
YXS   = p(3);
YBam  = p(4);
YAmX  = p(5);
YXO   = p(6);
KLa   = p(7);

% Kinetik und deren Ableitung nach cGlc
mu       = mumax * cGlc / (KS + cGlc);
dmudcGlc = mumax * KS / (KS + cGlc)^2;

% Spaltenweise Definition der Matrix (5x5)
dfdx = [
    mu,          dmudcGlc * cX,          0, 0, 0;          % d(dxdt1)/dx
    -mu / YXS,  -1/YXS * dmudcGlc * cX,  0, 0, 0;          % d(dxdt2)/dx
    -YAmX * mu, -YAmX * dmudcGlc * cX,   0, 0, 0;          % d(dxdt3)/dx
    YBam*YAmX*mu, YBam*YAmX*dmudcGlc*cX, 0, 0, 0;          % d(dxdt4)/dx
    -H * mu / YXO, -H * (1/YXO) * dmudcGlc * cX, 0, 0, -KLa% d(dxdt5)/dx
];
end