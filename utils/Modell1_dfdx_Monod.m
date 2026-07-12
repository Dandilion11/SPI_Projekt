function dfdx = Modell1_dfdx(x, p)
% Jacobimatrix df/dx fuer Modell1 (Monod, ohne O2)
% Zustände: x = [cX, cGlc], Parameter: p = [mumax, KS, YXS]

cX   = x(1);
cGlc = x(2);
mumax = p(1);
KS    = p(2);
YXS   = p(3);

mu        = mumax * cGlc / (KS + cGlc);
dmudcGlc  = mumax * KS  / (KS + cGlc)^2;

dfdx = [
           mu,          dmudcGlc * cX;
    -mu / YXS,  -1/YXS * dmudcGlc * cX
];
end