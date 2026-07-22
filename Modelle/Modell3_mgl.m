function y = Modell3_mgl(x)
% Messgleichung für Modell 3 (Krämer 2017, Gl. 20+22) ohne Ethanol
% ci = mi / V  für alle Massen-Zustände

V = x(1);

y = [x(2)/V;   % cX   [g/L]
     x(3)/V;   % cGlc [g/L]
     x(4)/V;   % cNH4 [g/L]
     x(5)/V;   % cPO4 [g/L]
     x(6);    % mB   [mL]
     x(7);];   % DOT [%]

end