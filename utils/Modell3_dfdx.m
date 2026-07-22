function dfdx = Modell3_dfdx(x, p, DOTstern)
% Jacobimatrix der Zustands-DGL nach den Zustaenden (Modell3 woEtOH)
x = max(x, 0);

V = x(1); mX = x(2); mGlc = x(3);

mumax = p(1); KS = p(2); YXS = p(3); YAmX = p(4);
YPhX  = p(5); YB_Am = p(6); KLa = p(7); YXO = p(8);

cO2_sat = 7.14e-3;
H = 100/cO2_sat;

cGlc = mGlc / V;
rX   = mumax * cGlc / (KS + cGlc);

% Ableitungen von rX
drdcGlc = mumax * KS / (KS + cGlc)^2;
drdG    = drdcGlc * (1/V);          % drX/dmGlc
drdV    = drdcGlc * (-mGlc/V^2);    % drX/dV

dfdx = zeros(7,7);
%        V                       mX            mGlc
% f1 = Zufluesse (nur Eingaenge) -> alle 0

% f2 = rX*mX
dfdx(2,1) = mX*drdV;   dfdx(2,2) = rX;             dfdx(2,3) = mX*drdG;

% f3 = -1/YXS*rX*mX
dfdx(3,1) = -1/YXS*mX*drdV;  dfdx(3,2) = -1/YXS*rX;  dfdx(3,3) = -1/YXS*mX*drdG;

% f4 = -YAmX*rX*mX
dfdx(4,1) = -YAmX*mX*drdV;   dfdx(4,2) = -YAmX*rX;   dfdx(4,3) = -YAmX*mX*drdG;

% f5 = -YPhX*rX*mX
dfdx(5,1) = -YPhX*mX*drdV;   dfdx(5,2) = -YPhX*rX;   dfdx(5,3) = -YPhX*mX*drdG;

% f6 = YB_Am*YAmX*rX*mX
dfdx(6,1) = YB_Am*YAmX*mX*drdV;  dfdx(6,2) = YB_Am*YAmX*rX;  dfdx(6,3) = YB_Am*YAmX*mX*drdG;

% f7 = KLa*(cO2stern-DOT) - H/YXO*rX*mX/V
dfdx(7,1) = -H/YXO * mX * (drdV/V - rX/V^2);
dfdx(7,2) = -H/YXO * rX/V;
dfdx(7,3) = -H/YXO * drdG * mX/V;
dfdx(7,7) = -KLa;
end