% gen_jacobians.m
%
% Erzeugt die symbolischen Ableitungen von Modell 3 als MATLAB-Funktionen:
%   Modell3_dfdx_10p.m   df/dx   (7 x 7)
%   Modell3_dfdp_10p.m   df/dp   (7 x 10)
%
% Beide werden fuer die Sensitivitaetsgleichungen gebraucht
% (Modell3_woEtOH_XP -> FIM):
%   dXP/dt = df/dx * XP + df/dp
%
% WICHTIG: f muss hier exakt dem RHS in Modell3_woEtOH_10p.m entsprechen.
% Nach jeder Modelaenderung dieses Skript neu laufen lassen und danach
% check_jacobians.m ausfuehren -- sonst rechnet die FIM mit einem anderen
% Modell als der Fit.
%
% Die Dateinamen muessen zu den Funktionsnamen passen: MATLAB waehlt beim
% Aufruf die DATEI. Ein Modell3_dfdp_10p.m, das intern Modell3_dfdp heisst,
% laesst sich zwar aufrufen, aber ein altes Modell3_dfdp.m auf dem Pfad
% wuerde still bevorzugt.

clear; clc;

syms V mX mGlc mAm mPh mB DOT real
syms mumax KS YXS YAmX YPhX YB_Am KLa YXO KAm KPh real
syms uGlc cGlc_in uAm cAm_in uPh cPh_in uBase uAcid DOTstern real

x = [V; mX; mGlc; mAm; mPh; mB; DOT];
p = [mumax; KS; YXS; YAmX; YPhX; YB_Am; KLa; YXO; KAm; KPh];

cO2_sat = 7.14e-3;   H = 100/cO2_sat;

% Dreifache Monod-Kinetik, ohne max(.,0): f muss differenzierbar sein.
cGlc = mGlc/V;   cAm = mAm/V;   cPh = mPh/V;
rX = mumax * (cGlc/(KS+cGlc)) * (cAm/(KAm+cAm)) * (cPh/(KPh+cPh));

f = [ uGlc + uAm + uPh + uBase + uAcid;                      % V
      rX*mX;                                                 % mX
      -1/YXS*rX*mX + cGlc_in*uGlc;                           % mGlc
      -YAmX*rX*mX  + cAm_in*uAm;                             % mAm
      -YPhX*rX*mX  + cPh_in*uPh;                             % mPh
      YB_Am*YAmX*rX*mX;                                      % mB
      KLa*(DOTstern-DOT) - H*(1/YXO)*rX*(mX/V) ];            % DOT

matlabFunction(jacobian(f,x), 'File','Modell3_dfdx_10p', ...
               'Vars',{x,p,DOTstern}, 'Optimize',true);
matlabFunction(jacobian(f,p), 'File','Modell3_dfdp_10p', ...
               'Vars',{x,p,DOTstern}, 'Optimize',true);

fprintf('Erzeugt: Modell3_dfdx_10p.m und Modell3_dfdp_10p.m\n');
fprintf('Jetzt check_jacobians.m ausfuehren.\n');