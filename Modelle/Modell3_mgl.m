function y = Modell3_mgl(x)
% Messgleichung für Modell 3 (Krämer 2017, Gl. 20+22)
% ci = mi / V  für alle Massen-Zustände

V = x(7);

y = [x(1)/V;   % cX   [g/L]
     x(2)/V;   % cGlc [g/L]
     x(3)/V;   % cNH4 [g/L]
     x(4)/V;   % cPO4 [g/L]
     x(5)/V;   % cEt  [g/L]
     x(6);    % mB   [mL]
     x(7);];   % Volumen

end