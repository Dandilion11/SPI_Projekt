function dxdt = Modell3(t, x, u, p)
% Modell 2c: Fed-Batch mit Glucose-, Ammonium- und Phosphat-Feed
% sowie Basezugabe. Zustände als Massen + Volumen.

%% Zustände
mX   = x(1);   % Biomasse
mS   = x(2);   % Glucose
mNH4 = x(3);   % Ammonium
mPO4 = x(4);   % Phosphat
mB   = x(5);   % akkumulierte Base
V    = x(6);   % Volumen

% keine negativen Zustände
x(x < 0) = 0;

%% Parameter
mumax  = p(1); % max spezifische Wachstumsrate (maximum specific growth rate (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
KS     = p(2); % Halbsättigungskonstante(substrate constant (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
YXS    = p(3); % Ertragskoeffizient (biomass yield (Doran, P. M. (2013). Bioprocess Engineering Principles, Gl. 4.11))
YNH4X  = p(4); 
YPO4X  = p(5);
YB_NH4 = p(6);
cS_F   = p(7);
cNH4_F = p(8);
cPO4_F = p(9);

%% Stellgrößen aus u (wie in Mess.u)
tu  = u(1,:);

idx = find(t >= tu, 1, 'last');  % aktueller Index

qC  = u(2, idx);   % C-Feedrate [L/h]
qAm = u(4, idx);   % Ammonium-Feedrate [L/h]
qPh = u(6, idx);   % Phosphat-Feedrate [L/h]
qB  = u(10, idx);  % Basezugabe [L/h]

cS = mS / V;   % Konzentration im Reaktor[g/L]

%% Wachstumsrate (Monod)
mu = mumax * cS / (KS + cS);

%% Differentialgleichungen
dxdt = zeros(6,1);

% Biomasse
dxdt(1) = mu * mX;

% Glucose
dxdt(2) = - (1 / YXS) * mu * mX + cS_F   * qC;

% Ammonium
dxdt(3) = - (1 / YNH4X) * mu * mX + cNH4_F * qAm;

% Phosphat
dxdt(4) = - (1 / YPO4X) * mu * mX + cPO4_F * qPh;

% Base (Feed + proportional zur Ammoniumassimilation)
dxdt(5) = qB + YB_NH4 * ( -dxdt(3) );

% Volumen (Summe der Volumenströme)
dxdt(6) = qC + qAm + qPh + qB;
