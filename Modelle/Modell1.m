% DGL einfachem Modell (Aufgabe 2a) ohne Inputs.
% kinetic 1 für Moser/Blachmann, Kinetic 2 für exponential und Kinetic 3
% für Monod

function dxdt = Modell1(t, x, p, kinetic)

x = max(x, 0);

%% Zustände
cX   = x(1);  % Biomasse
cGlc = x(2);  % Glucose
cO2  = x(3);  % Sauerstoff

%% Parameter
mumax = p(1); % max spezifische Wachstumsrate (maximum specific growth rate (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
KS    = p(2); % Halbsättigungskonstante(substrate constant (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
YXS   = p(3); % Ertragskoeffizient Biomasse/Glucose (biomass yield (Doran, P. M. (2013). Bioprocess Engineering Principles, Gl. 4.11))  
YXO   = p(4);  % Ertragskoeffizient Biomasse/O2
KLa   = p(5);  % Aus Enfors Gl. 6.6 Volumetrischer Sauerstofftransferkoeffizient [1/h]
cO2stern = p(6);  % Enfors Gl. 6.6 -> cO2stern(C*): dissolved oxygen concentration in equilibrium with the gas phase -> Maximal lösliche O2-Konzentration unter den gegebenen Bedingungen

%% Wachstumsrate
if kinetic == 1
    n  = p(7);
    mu = mumax * (cGlc^n / (cGlc^n + KS));
elseif kinetic == 2
    mu = mumax * (1 - exp(-cGlc / KS));
else
    mu = mumax * cGlc / (KS + cGlc);
end

%% Differentialgleichungen
dxdt = zeros(3, 1);

dxdt(1) =  mu * cX;                                 % Biomasse
dxdt(2) = -1/YXS * mu * cX;                         % Glucose
dxdt(3) =  KLa * (cO2stern-cO2) - 1/YXO * mu * cX;  % O2-Eintrag - O2-Verbrauch (Enfors Gl. 6.15) -> O2 Eintrag durch den Rührer, O2 Verbrauch durch das Wachstum der Biomasse 
 
% in Enfors: d(DOT)/dt = KLa * (DOT* - DOT) - r0*H  (Gl. 6.15)
% -> DOT = C * 100 kH / pO2,cal                     (Gl. 6.14 (umgestellt))
% -> H = 100 kH / PO2,cal              (nach Gl. 6.15 im Text gegeben)
% -> DOT = C * H = cO2 * H    (C zu cO2 in unserem Beispiel definiert)
% -> cO2 = DOT/H
% -> d(cO2)/dt = KLa * (cO2stern - cO2) - r0 (cO2 mit DOT/H ersetzt)
% -> r0 = q0 * cX = mu/YXO * cX (Gl. 6.1)
% ---> d(cO2)/dt = KLa(cO2stern-cO2) - 1/YXO * mu * cX




end