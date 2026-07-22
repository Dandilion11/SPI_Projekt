function dhdx = Modell3_dmgldx(x)
% Jacobimatrix der Messgleichungen h nach den Zustaenden x.
% h = [cX; cGlc; cAm; cPh; mB; DOT]   (6 Messgroessen)
% x = [V; mX; mGlc; mAm; mPh; mB; DOT] (7 Zustaende)
% dhdx : 6 x 7

V   = x(1);
mX  = x(2);
mGlc= x(3);
mAm = x(4);
mPh = x(5);

dhdx = zeros(6,7);

% h1 = cX = mX/V
dhdx(1,1) = -mX / V^2;    % d/dV
dhdx(1,2) =  1  / V;      % d/dmX

% h2 = cGlc = mGlc/V
dhdx(2,1) = -mGlc / V^2;
dhdx(2,3) =  1    / V;

% h3 = cAm = mAm/V
dhdx(3,1) = -mAm / V^2;
dhdx(3,4) =  1   / V;

% h4 = cPh = mPh/V
dhdx(4,1) = -mPh / V^2;
dhdx(4,5) =  1   / V;

% h5 = mB (direkt gemessen)
dhdx(5,6) = 1;

% h6 = DOT (direkt gemessen)
dhdx(6,7) = 1;
end