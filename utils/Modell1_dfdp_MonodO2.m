function dfdp = Modell1_dfdp_MonodO2(x, p)
% Jacobimatrix df/dp für Modell1 (Monod) mit Sauerstoff
% Ableitungen erstmal nur für Monod kinetik berechnet

% Prinzip Sensitivitäten für nichtlineare Modelle: df/dx * inv(C) + df/dp
% hier wird df/dp berechnet
% im Gegensatz zu Übung 5 wird hier eine 3x5 matrix rauskommen, weil ohne
% Volumen gerechnet wird, weil nur von Konzentrationen ausgegangen wird

cO2stern = 100;  % Enfors Gl. 6.6 -> cO2stern(C*): dissolved oxygen concentration in equilibrium with the gas phase -> Maximal lösliche O2-Konzentration unter den gegebenen Bedingungen
cO2_sat = 7.14e-3;   % Laut Enfors für das Modell angenommen (g/L)
H       = 100/cO2_sat; 

cX   = x(1);
cGlc = x(2);
cO2 = x(3);

mumax = p(1);
KS    = p(2);
YXS   = p(3);  
YXO = p(4);


mu      = mumax * cGlc / (KS + cGlc);
dmumax  = cGlc / (KS + cGlc);              % dmu/d(mumax)
dmudKS  = -mumax * cGlc / (KS + cGlc)^2;  % dmu/d(KS)

% Zeilen x1 -> x3
% Spalten: d/dmumax     d/dKS             d/dYXS       d/dYXO         d/KLa
dfdp = [
       dmumax*cX,       dmudKS*cX,        0,           0,             0;
      -dmumax*cX/YXS,   -dmudKS*cX/YXS,   mu*cX/YXS^2, 0,             0;
      -H/YXO*dmumax*cX, -H/YXO*dmudKS*cX, 0,           H*mu*cX/YXO^2, cO2stern - cO2
    ];
end