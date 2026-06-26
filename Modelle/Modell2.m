function dxdt = Modell2(t,x,p)
% Erweiterung des Modells 1 um eine base consumption:
% für ph regulation eines säure-produzierenden Prozesses
%% Grenzen der Zustände

x = max(x,0);

%% Zustände

cX   = x(1);
cS   = x(2);
cNH4 = x(3); %Ammonium
cB   = x(4); %akkumulierte Base
%% Parameter

mumax = p(1); % max spezifische Wachstumsrate (maximum specific growth rate (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
KS    = p(2); % Halbsättigungskonstante(substrate constant (Doran, P. M. (2013). Bioprocess Engineering Principles. S. 279))
YXS   = p(3); % Ertragskoeffizient (biomass yield (Doran, P. M. (2013). Bioprocess Engineering Principles, Gl. 4.11))  
YBam  = p(4);  % yield Base/Ammonium 
YAmX  = p(5);  % yield Ammonium/Biomasse

%% Wachstumsrate

mu = mumax * cS / (KS + cS);

%% DGL

dxdt = zeros(4,1);

dxdt(1) = mu * cX; % Biomasse
dxdt(2) = - (1 / YXS) * mu * cX; % Glucose
dxdt(3) = - YAmX * mu * cX; % Ammonium
dxdt(4) = YBam * YAmX * mu * cX; %Base

end