function dxdt = Modell3_YBA(t, x, u, p, DOTstern)
% Modell 3: Fed-Batch MIT Ethanol (Krämer & King 2017, Gl. 7–15)
% Sauerstoff nicht als Zustand, weil der Konstant auf >=50% gehalten wird,
% durch den Rührer
%

x = max(x, 0);

% Nach Input struct     % Zustände entsprechend des
% aus Preprocessing     % Mess.Messdaten-Structs                 
V    = x(1);            % V    = x(1);
mX   = x(2);            % mX   = x(3);
mGlc = x(3);            % mGlc = x(4);
mAm  = x(4);            % mAm  = x(5);
mPh  = x(5);            % mPh  = x(6);
mB   = x(6);            % mB   = x(7);
DOT  = x(7);            % DOT  = x(8);
mEt  = x(8);            % mEt  = x(9);

cO2stern = DOTstern; %Laut Nachricht von Terrance max(DOT) statt 100 aus Enfors -> Enfors Gl. 6.6 -> cO2stern(C*): dissolved oxygen concentration in equilibrium with the gas phase -> Maximal lösliche O2-Konzentration unter den gegebenen Bedingungen
cO2_sat = 7.14e-3;   % Laut Enfors für das Modell angenommen (g/L)
H       = 100/cO2_sat; 
%% Parameter
mumax     = p(1);
KS        = p(2);
YXS       = p(3);
YPhX      = p(4);
YBA       = p(5);
mumax_EtP = p(6);
mumax_EtX = p(7);
YGlc_Et   = p(8);
YEt_X     = p(9);
KEt       = p(10);
KGlc_Et   = p(11);
KLa       = p(12);
YXO       = p(13);

%% Konzentrationen
cGlc = mGlc / V;
cEt  = mEt  / V;
cAm = mAm / V;
cPh = mPh / V;

%% Reaktionsraten (Krämer 2017, Gl. 16–18)

% rX: Biomassewachstum auf Glucose (Gl. 16)                   % Monod (Standard)
 rX = mumax * cGlc / (KS + cGlc);

% rEt_P: Ethanolproduktion (Crabtree-Effekt, Gl. 17): Aktiv bei hoher Glukose → Overflow-Metabolismus
rEt_P = mumax_EtP * cGlc / (cGlc + KS);

% rEt_X: Ethanolverbrauch (Gl. 18): Aktiv bei niedriger Glukose + vorhandenem Ethanol
rEt_X = mumax_EtX * (cEt / (cEt + KEt)) * (KGlc_Et / (cGlc + KGlc_Et));

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
dxdt = zeros(8, 1);
dxdt(1) = uGlc + uAm + uPh + uBase + uAcid - uout;                      % V
dxdt(2) = (rX + rEt_X) * mX;                                            % mX
dxdt(3) = (-1/YXS*rX - YGlc_Et*rEt_P)*mX + cGlc_in*uGlc - uout*cGlc;    % mGlc
dxdt(4) = -YAmX*(rX+rEt_X)*mX + cAm_in*uAm - uout*cAm;                  % mAm
dxdt(5) = -YPhX*(rX+rEt_X)*mX + cPh_in*uPh - uout*cPh;                  % mPh
dxdt(6) = YBA *(rX+rEt_X)*mX;                                           % mB YB_Am*YAmX -> YBA
dxdt(7) = KLa*(cO2stern-DOT) - H * (1/YXO)*(rX+rEt_X)*(mX/V);           % DOT
dxdt(8) = (rEt_P - YEt_X*rEt_X)*mX - uout*cEt;                          % mEt 

end

