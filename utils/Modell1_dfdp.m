function dfdp = Modell1_dfdp(x, p)
% Jacobimatrix df/dp für Modell1 (Monod, ohne Sauerstoff)

cX   = x(1);
cGlc = x(2);
mumax = p(1);
KS    = p(2);
YXS   = p(3);

mu      = mumax * cGlc / (KS + cGlc);
dmumax  = cGlc / (KS + cGlc);              % dmu/d(mumax)
dmudKS  = -mumax * cGlc / (KS + cGlc)^2;  % dmu/d(KS)

%         d/dmumax          d/dKS           d/dYXS
dfdp = [
         dmumax*cX,       dmudKS*cX,             0;
    -dmumax*cX/YXS,  -dmudKS*cX/YXS,  mu*cX/YXS^2
];
end