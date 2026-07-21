function dfdx = Modell1_dfdx_MonodO2(x, p)
% Jacobimatrix df/dx für Modell1 (Monod, mit O2)
% Prinzip Sensitivitäten für nichtlineare Modelle: df/dx * inv(C) + df/dp
% im Gegensatz zu Übung 5 wird hier eine 3x5 matrix rauskommen, weil ohne
% Volumen gerechnet wird, weil nur von Konzentrationen ausgegangen wird

cO2_sat = 7.14e-3;      % g/L
H       = 100 / cO2_sat;

cX   = x(1);
cGlc = x(2);
mumax = p(1);
KS    = p(2);
YXS   = p(3);
YXO = p(4);
KLa = p(5);

mu        = mumax * cGlc / (KS + cGlc);
dmudcGlc  = mumax * KS  / (KS + cGlc)^2;

dfdx = [
                mu,              dmudcGlc * cX,          0;
         -mu / YXS,      -1/YXS * dmudcGlc * cX,          0;
        -H/YXO*mu,       -H/YXO * dmudcGlc * cX,       -KLa
    ];
end