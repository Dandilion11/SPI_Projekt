function dxdt = Modell3(t, x, u, p, kinetic)
% Modell 3: Fed-Batch MIT Ethanol (Krämer & King 2017, Gl. 7–15)
%
% Zustände (Massen [g] + Volumen [L]):
%   x(1) = mX    Biomasse
%   x(2) = mGlc  Glucose
%   x(3) = mNH4  Ammonium
%   x(4) = mPO4  Phosphat
%   x(5) = mEt   Ethanol    
%   x(6) = mB    kumulierte Base [mL]
%   x(7) = V     Volumen [L]
%
% Parameter p:
%   p(1)  mumax      max. Wachstumsrate auf Glucose [1/h]
%   p(2)  KS         Halbsättigungskonstante Glucose [g/L]
%   p(3)  YXS        Yield Biomasse/Glucose [g_X/g_Glc]
%   p(4)  YNH4X      Yield NH4/Biomasse [g_NH4/g_X]
%   p(5)  YPO4X      Yield PO4/Biomasse [g_PO4/g_X]
%   p(6)  YB_NH4     Yield Base/NH4 [mL/g_NH4]     (Krämer Gl. 13)
%   p(7)  mumax_EtP  max. Ethanolproduktionsrate [1/h]
%   p(8)  mumax_EtX  max. Wachstumsrate auf Ethanol [1/h]
%   p(9)  YGlc_Et    Yield Glucose/EtOH-Produktion [g_Glc/g_EtOH]
%   p(10) YEt_X      Yield Ethanol/Biomasse (Konsum) [g_EtOH/g_X]
%   p(11) KEt        Sättigungskonstante Ethanol [g/L]
%   p(12) KGlc_Et    Inhibierung Glc → EtOH-Konsum [g/L]
%   p(13) n          Moser-Exponent (nur wenn kinetic == 1)

x = max(x, 0);

%% Zustände
mX   = x(1);
mGlc = x(2);
mNH4 = x(3);
mPO4 = x(4);
mEt  = x(5);
V    = x(6);

%% Parameter
mumax     = p(1);
KS        = p(2);
YXS       = p(3);
YNH4X     = p(4);
YPO4X     = p(5);
YB_NH4    = p(6);
mumax_EtP = p(7);
mumax_EtX = p(8);
YGlc_Et   = p(9);
YEt_X     = p(10);
KEt       = p(11);
KGlc_Et   = p(12);

%% Konzentrationen
cGlc = mGlc / V;
cEt  = mEt  / V;
cNH4 = mNH4 / V;
cPO4 = mPO4 / V;

%% Reaktionsraten (Krämer 2017, Gl. 16–18)

% rX: Biomassewachstum auf Glucose (Gl. 16)
if kinetic == 1          % Moser
    n  = p(13);
    rX = mumax * (cGlc^n / (cGlc^n + KS));
elseif kinetic == 2      % Tessier
    rX = mumax * (1 - exp(-cGlc / KS));
else                     % Monod (Standard)
    rX = mumax * cGlc / (KS + cGlc);
end

% rEt_P: Ethanolproduktion (Crabtree-Effekt, Gl. 17)
% Aktiv bei hoher Glukose → Overflow-Metabolismus
rEt_P = mumax_EtP * cGlc / (cGlc + KS);

% rEt_X: Ethanolverbrauch (Gl. 18)
% Aktiv bei niedriger Glukose + vorhandenem Ethanol
rEt_X = mumax_EtX * (cEt / (cEt + KEt)) * (KGlc_Et / (cGlc + KGlc_Et));

%% Stellgrößen aus u-Matrix (Krämer 2016/2017, plotter_u.m)
tu  = u(1,:);
idx = find(t >= tu, 1, 'last');

qAm   = u(2,  idx);   % Ammonium-Feedrate [L/h]
cNH4F = u(3,  idx);   % NH4-Konz. im Feed [g/L]
qPh   = u(4,  idx);   % Phosphat-Feedrate [L/h]
cPO4F = u(5,  idx);   % PO4-Konz. im Feed [g/L]
qGlc  = u(6,  idx);   % Glucose-Feedrate [L/h]
cGlcF = u(7,  idx);   % Glc-Konz. im Feed [g/L]  (≈450 g/L)
qAcid = u(9,  idx);   % Säure-Feedrate [L/h] Säure Feedrate = 0, also nicht relevant, aber trotzdem mitbetrachtet, der vollständigkeitshalber
qBase = u(10, idx);   % Base-Feedrate [L/h]

%% Probenahme-Abfluss
% Probenahmen diskret (~1–2 mL) → uout = 0 in der DGL.
% Diskrete Zustandskorrektur für Probenentnahme in Simulation
uout = 0;

%% Massenbilanzen (Krämer & King 2017, Gl. 7–15)
dxdt = zeros(7, 1);

% Biomasse (Gl. 7): Wachstum auf Glc + Wachstum auf EtOH
dxdt(1) = (rX + rEt_X) * mX;

% Glucose (Gl. 10): Verbrauch durch Wachstum + EtOH-Produktion + Feed
dxdt(2) = (-1/YXS * rX - YGlc_Et * rEt_P) * mX ...
          + cGlcF * qGlc - uout * cGlc;

% Ammonium (Gl. 8): Verbrauch gesamt + Feed
dxdt(3) = -YNH4X * (rX + rEt_X) * mX ...
          + cNH4F * qAm - uout * cNH4; 
% Phosphat (Gl. 9): Verbrauch gesamt + Feed
dxdt(4) = -YPO4X * (rX + rEt_X) * mX ...
          + cPO4F * qPh - uout * cPO4; 

% Ethanol (Gl. 11): Produktion – Verbrauch
dxdt(5) = (rEt_P - YEt_X * rEt_X) * mX - uout * cEt;    

% Kumulierte Base (Gl. 13)
dxdt(6) = YB_NH4 * YNH4X * (rX + rEt_X) * mX;

% Volumen (Gl. 15)
dxdt(7) = qGlc + qAm + qPh + qBase + qAcid - uout;   

end

