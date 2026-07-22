function dfdp = Modell3_dfdp(x, p, DOTstern)
% Ableitung der Zustands-DGL nach den Parametern (Modell3 woEtOH)
x = max(x, 0);

V = x(1); mX = x(2); mGlc = x(3); DOT = x(7);

mumax = p(1); KS = p(2); YXS = p(3); YAmX = p(4);
YPhX  = p(5); YB_Am = p(6); KLa = p(7); YXO = p(8);

cO2_sat  = 7.14e-3;
H        = 100/cO2_sat;
cO2stern = DOTstern;

cGlc = mGlc / V;
rX   = mumax * cGlc / (KS + cGlc);

% Ableitungen von rX nach Parametern
drdmu = cGlc / (KS + cGlc);
drdKS = -mumax * cGlc / (KS + cGlc)^2;

dfdp = zeros(7,8);
% Spalten: mumax KS YXS YAmX YPhX YB_Am KLa YXO

% f2 = rX*mX
dfdp(2,1) = mX*drdmu;      dfdp(2,2) = mX*drdKS;

% f3 = -1/YXS*rX*mX
dfdp(3,1) = -1/YXS*mX*drdmu;   dfdp(3,2) = -1/YXS*mX*drdKS;   dfdp(3,3) = 1/YXS^2*rX*mX;

% f4 = -YAmX*rX*mX
dfdp(4,1) = -YAmX*mX*drdmu;    dfdp(4,2) = -YAmX*mX*drdKS;    dfdp(4,4) = -rX*mX;

% f5 = -YPhX*rX*mX
dfdp(5,1) = -YPhX*mX*drdmu;    dfdp(5,2) = -YPhX*mX*drdKS;    dfdp(5,5) = -rX*mX;

% f6 = YB_Am*YAmX*rX*mX
dfdp(6,1) = YB_Am*YAmX*mX*drdmu;  dfdp(6,2) = YB_Am*YAmX*mX*drdKS;
dfdp(6,4) = YB_Am*rX*mX;          dfdp(6,6) = YAmX*rX*mX;

% f7 = KLa*(cO2stern-DOT) - H/YXO*rX*mX/V
dfdp(7,1) = -H/YXO*drdmu*mX/V;
dfdp(7,2) = -H/YXO*drdKS*mX/V;
dfdp(7,7) = (cO2stern - DOT);
dfdp(7,8) =  H/YXO^2 * rX * mX/V;
end