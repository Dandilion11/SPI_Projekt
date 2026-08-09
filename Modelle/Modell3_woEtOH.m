function dxdt = Modell3_woEtOH(t, x, u, p, DOTstern)
% Modell 3: Fed-Batch ohne (Krämer & King 2017, Gl. 7–15)

%Ethanol ist in den Messwerten immer 0 und wurde nicht über die ganze
%Messdauer aufgezeichnet deshalb die Version ohne EtOH

%x = max(x, 0);

% Nach Input struct     % Zustände entsprechend des
% aus Preprocessing     % Mess.Messdaten-Structs                 
V    = x(1);            % V    = x(1);
mX   = x(2);            % mX   = x(3);
mGlc = x(3);            % mGlc = x(4);
mAm  = x(4);            % mAm  = x(5);
mPh  = x(5);            % mPh  = x(6);
mB   = x(6);            % mB   = x(7);
DOT  = x(7);            % DOT  = x(8);

cO2stern = DOTstern; %Laut Nachricht von Terrance max(DOT) statt 100 aus Enfors -> Enfors Gl. 6.6 -> cO2stern(C*): dissolved oxygen concentration in equilibrium with the gas phase -> Maximal lösliche O2-Konzentration unter den gegebenen Bedingungen
cO2_sat = 7.14e-3;   % Laut Enfors für das Modell angenommen (g/L)
H       = 100/cO2_sat; 
%% Parameter
mumax     = p(1);
KS        = p(2);
YXS       = p(3);
YAmX      = p(4);
YPhX      = p(5);
YB_Am     = p(6);
KLa       = p(7);
YXO       = p(8);

%% Konzentrationen
cGlc = mGlc / V;
cAm = mAm / V;
cPh = mPh / V;

%% Reaktionsraten (Krämer 2017, Gl. 16–18)

% rX: Biomassewachstum auf Glucose (Gl. 16)                   % Monod (Standard)
 rX = mumax * (cGlc/(KS + cGlc));  % [1/h]

%% Stellgrößen aus u-Matrix (Krämer 2016/2017 und laut Terrance Beschreibung der Eingänge)
tu  = u(1,:);
idx = find(t >= tu, 1, 'last');
if isempty(idx)
    idx = 1;                 % vor erstem Zeitpunkt: ersten Wert halten
end
idx = min(idx, size(u,2));

uAm     = u(2,  idx);   % Ammonium-Feedrate [L/h] 
cAm_in  = u(3,  idx);  % Ammonium-Eingangskonzentration [g/L] (30 g/L)

uPh     = u(4,  idx);   % Phosphat-Feedrate [L/h] 
cPh_in  = u(5,  idx);   % Phosphat-Eingangskonzentration [g/L] (24 g/L)

uGlc    = u(6,  idx);   % Glucose-Feedrate [L/h] 
cGlc_in = u(7,  idx);   % Glc-Konz. im Feed [g/L] (450 g/L)

uAcid   = u(9,  idx);   % Säure-Feedrate [L/h]
uBase   = u(10, idx);   % Base-Feedrate [L/h]

%% Probenahme-Abfluss
% Diskrete Zustandskorrektur für Probenentnahme in Simulation
uout = 0;

%% Massenbilanzen (Krämer & King 2017, Gl. 7–15)
dxdt = zeros(7, 1);
dxdt(1) = uGlc + uAm + uPh + uBase + uAcid - uout;              % V
dxdt(2) = (rX) * mX;                                            % mX
dxdt(3) = (-1/YXS*rX)*mX + cGlc_in*uGlc - uout*cGlc;            % mGlc
dxdt(4) = -YAmX*rX*mX + cAm_in*uAm - uout*cAm;                  % mAm
dxdt(5) = -YPhX*rX*mX + cPh_in*uPh - uout*cPh;                  % mPh
dxdt(6) = YB_Am*YAmX*rX*mX;                                     % mB
dxdt(7) = KLa*(cO2stern-DOT) - H * (1/YXO)*(rX)*(mX/V);         % DOT

end
