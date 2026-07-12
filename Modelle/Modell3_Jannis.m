function dxdt = Modell2(t, x, p, u, DOTstern)
% Modell2 angepasst an das vollständige Modell:
%
% x = [V; m_X; m_glc; DOT; m_am; m_ph; m_b]
% p = [mu_glc_max_X; K_glc_X; Y_glc_X; Kla; Y_X_O2; Y_am_X; Y_ph_X; Y_B_am]
% u = [u_glc_in; c_glc_in; u_am_in; c_am_in; u_ph_in; c_ph_in; u_B_in; c_B_in]


u = u(:);
x = max(x(:), 0);


%% Zustände
V     = x(1);
m_X   = x(2);  % Biomasse
m_glc = x(3);  % Glucose
DOT   = x(4);  % Dissolved Oxygen Tension
m_am  = x(5);  % Ammonium
m_ph  = x(6);  % Phosphat
m_b   = x(7);  % Base

%% Konzentrationen
c_X   = m_X / V;
c_glc = m_glc / V;

%% Parameter
mu_glc_max_X = p(1);
K_glc_X      = p(2);
Y_glc_X      = p(3);
Kla          = p(4);
Y_X_O2       = p(5);
Y_am_X       = p(6);
Y_ph_X       = p(7);
Y_B_am       = p(8);

%% Inputs
u_glc_in = u(1);
c_glc_in = u(2);

u_am_in  = u(3);
c_am_in  = u(4);

u_ph_in  = u(5);
c_ph_in  = u(6);

u_B_in   = u(7);
c_B_in   = u(8);

%% Sauerstoff-Konstanten
DOTstar = DOTstern;

cO2_sat = 7.14e-3;
H = 100 / cO2_sat;

%% Reaktionsrate / Monod-Kinetik
r_X = (mu_glc_max_X * c_glc) / (K_glc_X + c_glc);

%% DGL-System
dxdt = zeros(7, 1);

dxdt(1) = u_glc_in + u_am_in + u_ph_in + u_B_in;

dxdt(2) = r_X * m_X;

dxdt(3) = -Y_glc_X * r_X * m_X ...
          + u_glc_in * c_glc_in;

dxdt(4) = Kla * H * (DOTstar - DOT) ...
          - c_X * r_X * Y_X_O2;

dxdt(5) = -Y_am_X * r_X * m_X ...
          + u_am_in * c_am_in;

dxdt(6) = -Y_ph_X * r_X * m_X ...
          + u_ph_in * c_ph_in;

dxdt(7) = Y_B_am * Y_am_X * r_X * m_X ...
          + u_B_in * c_B_in;

end