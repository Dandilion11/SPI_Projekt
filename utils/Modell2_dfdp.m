function dfdp = Modell2_dfdp(x, p)
% Jacobimatrix df/dp für Modell2 (Monod-Kinetik)
% Zustände: x = [cX, cGlc, cAm, cBase, cO2]
% Parameter: p = [mumax, KS, YXS, YBam, YAmX, YXO, KLa, cO2stern]

cO2stern = 100;  % Enfors Gl. 6.6 -> cO2stern(C*): dissolved oxygen concentration in equilibrium with the gas phase -> Maximal lösliche O2-Konzentration unter den gegebenen Bedingungen
cO2_sat = 7.14e-3;   % Laut Enfors für das Modell angenommen (g/L)
H       = 100/cO2_sat; 

cX       = x(1);
cGlc     = x(2);
cO2      = x(5);

mumax    = p(1);
KS       = p(2);
YXS      = p(3);
YBam     = p(4);
YAmX     = p(5);
YXO      = p(6);
KLa      = p(7);

mu     = mumax * cGlc / (KS + cGlc);
dmu_dm = cGlc / (KS + cGlc);             % dmu/d(mumax)
dmu_dK = -mumax * cGlc / (KS + cGlc)^2;  % dmu/d(KS)

% Matrix-Dimension: 5 Zustände x 8 Parameter
dfdp = zeros(5, 8);

% Parameter 1: mumax
dfdp(1,1) = dmu_dm * cX;
dfdp(2,1) = -1/YXS * dmu_dm * cX;
dfdp(3,1) = -YAmX * dmu_dm * cX;
dfdp(4,1) = YBam * YAmX * dmu_dm * cX;
dfdp(5,1) = -1/YXO * dmu_dm * cX;

% Parameter 2: KS
dfdp(1,2) = dmu_dK * cX;
dfdp(2,2) = -1/YXS * dmu_dK * cX;
dfdp(3,2) = -YAmX * dmu_dK * cX;
dfdp(4,2) = YBam * YAmX * dmu_dK * cX;
dfdp(5,2) = H * -1/YXO * dmu_dK * cX;

% Parameter 3: YXS
dfdp(2,3) = (1 / YXS^2) * mu * cX;

% Parameter 4: YBam
dfdp(4,4) = YAmX * mu * cX;

% Parameter 5: YAmX
dfdp(3,5) = -mu * cX;
dfdp(4,5) = YBam * mu * cX;

% Parameter 6: YXO
dfdp(5,6) = (H / YXO^2) * mu * cX;

% Parameter 7: KLa
dfdp(5,7) = cO2stern - cO2;

end