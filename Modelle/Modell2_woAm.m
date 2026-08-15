function dxdt = Modell2_woAm(t,x,p, kinetic, DOTstern)
% Erweiterung des Modells 1 um eine Base-Consumption (pH-Regelung eines
% saeure-/ammonium-zehrenden Prozesses) sowie eine Sauerstoffbilanz.
%% Grenzen der Zustände

x = max(x,0);


%% Zustaende
cX    = x(1);  % Biomasse             [g/L] (CDW)
cGlc  = x(2);  % Glucose              [g/L]
cBase = x(3);  % zugegebene Base      [mol/L bzw. L-Aequivalent]
cO2   = x(4);  % gelöster Sauerstoff [g/L]   

% Enfors Gl. 6.6 -> cO2stern(C*): dissolved oxygen concentration in equilibrium with the gas phase -> Maximal lösliche O2-Konzentration unter den gegebenen Bedingungen
cO2_sat = 7.14e-3;   % Laut Enfors für das Modell angenommen (g/L)
H       = 100/cO2_sat; 
%% Parameter

mumax = p(1);       % max spezifische Wachstumsrate (maximum specific growth rate (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
KS    = p(2);       % Halbsättigungskonstante(substrate constant (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
YXS   = p(3);       % Ertragskoeffizient (biomass yield (Doran, P. M. (2013). Bioprocess Engineering Principles, Gl. 4.11))  
YBA   = p(4);       % yield Base/Ammonium 
YXO   = p(5);  % Ertragskoeffizient Biomasse/O2
KLa   = p(6);  % Aus Enfors Gl. 6.6 Volumetrischer Sauerstofftransferkoeffizient [1/h]

%% Wachstumsrate
if kinetic == 1
    n = p(8);
    mu = mumax * (cGlc^n / (cGlc^n + KS));
elseif kinetic == 2
    mu = mumax * (1 - exp(-cGlc / KS));
else
    mu = mumax * cGlc / (KS + cGlc);
end


%% DGL
dxdt = zeros(4, 1);
dxdt(1) =  mu * cX;                                 % Biomasse
dxdt(2) = -(1 / YXS) * mu * cX;                     % Glucose
dxdt(3) =  YBA * mu * cX;                   % Base      ->Messung YBA ist Produkt aus YBam und YAmX
dxdt(4) =  KLa * (DOTstern - cO2) - H * 1/YXO * mu * cX; % O2   ->Messung Basen Massen Bilanz

end