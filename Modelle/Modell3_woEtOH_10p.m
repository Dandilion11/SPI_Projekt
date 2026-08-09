function dxdt = Modell3_woEtOH_10p(t, x, u, p, DOTstern)
% Modell 3: Fed-Batch ohne Ethanol (Kraemer & King 2017, Gl. 7-15).
%
% Ethanol ist in den Messwerten immer 0 und wurde nicht ueber die ganze
% Messdauer aufgezeichnet -- deshalb die Version ohne EtOH.
%
% Zustaende:  x = [V; mX; mGlc; mAm; mPh; mB; DOT]
% Parameter:  p = [mumax; KS; YXS; YAmX; YPhX; YB_Am; KLa; YXO; KAm; KPh]
%
% KAm und KPh sind gegenueber der Urfassung ergaenzt: ohne sie laeuft rX
% weiter, wenn Ammonium oder Phosphat aufgebraucht sind, und die
% Massenbilanz zieht mAm/mPh unter null. Gefittet werden sie nicht.
%
% Keine max(.,0)-Klammerung der Zustaende: die Sensitivitaetsgleichungen
% (Modell3_dfdx_10p / _dfdp_10p) setzen ein differenzierbares f voraus.
% Positivitaet sichert stattdessen die NonNegative-Option des Integrators.

%% Zustaende
V    = x(1);
mX   = x(2);
mGlc = x(3);
mAm  = x(4);
mPh  = x(5);
DOT  = x(7);        % mB = x(6) geht in kein Ratengesetz ein

%% Konstanten
cO2stern = DOTstern;      % max(DOT) aus den Daten statt 100 (Enfors Gl. 6.6)
cO2_sat  = 7.14e-3;       % g/L, Annahme nach Enfors
H        = 100/cO2_sat;   % ~1.4006e4 % L/g

%% Parameter
mumax = p(1);   KS    = p(2);   YXS = p(3);   YAmX = p(4);   YPhX = p(5);
YB_Am = p(6);   KLa   = p(7);   YXO = p(8);   KAm  = p(9);   KPh  = p(10);

%% Wachstumsrate: dreifache Monod-Kinetik (Kraemer 2017, Gl. 16)
cGlc = mGlc / V;
cAm  = mAm  / V;
cPh  = mPh  / V;

rX = mumax * (cGlc/(KS+cGlc)) * (cAm/(KAm+cAm)) * (cPh/(KPh+cPh));   % 1/h

%% Stellgroessen: aktueller Abschnitt der u-Matrix (Nullter-Ordnung-Halt)
idx = find(t >= u(1,:), 1, 'last');
if isempty(idx), idx = 1; end          % vor dem ersten Zeitpunkt
idx = min(idx, size(u,2));

uAm     = u(2, idx);   cAm_in  = u(3, idx);   % Ammonium-Feed  [L/h], [g/L]
uPh     = u(4, idx);   cPh_in  = u(5, idx);   % Phosphat-Feed  [L/h], [g/L]
uGlc    = u(6, idx);   cGlc_in = u(7, idx);   % Glucose-Feed   [L/h], [g/L]
uAcid   = u(9, idx);                          % Saeure         [L/h]
uBase   = u(10,idx);                          % Base           [L/h]

%% Massenbilanzen (Kraemer & King 2017, Gl. 7-15)
% Die Probenahme wird nicht hier, sondern als Sprung zwischen den
% Integrationssegmenten behandelt (probe_m3 in sim_m3_sample_10p).
dxdt = zeros(7,1);
dxdt(1) = uGlc + uAm + uPh + uBase + uAcid;                % V
dxdt(2) =  rX * mX;                                        % mX
dxdt(3) = -rX * mX / YXS   + cGlc_in*uGlc;                 % mGlc
dxdt(4) = -rX * mX * YAmX  + cAm_in *uAm;                  % mAm
dxdt(5) = -rX * mX * YPhX  + cPh_in *uPh;                  % mPh
dxdt(6) =  rX * mX * YAmX * YB_Am;                         % mB
dxdt(7) = KLa*(cO2stern - DOT) - H * rX * (mX/V) / YXO;    % DOT
end