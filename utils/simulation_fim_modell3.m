%% FIM von Modell 3
function [FM, x_out, y_out] = simulation_fim_modell3(t, x0, u, p, invC, DOTstern)
% Simulation mit Sensitivitaeten und Berechnung der Fisher-Informationsmatrix
% invC : 6x6 (inverse Messkovarianz, 6 Messgroessen)
nx = 7; np = 8; ny = 6;

if nargin < 6 || isempty(DOTstern), DOTstern = 100; end

% Anfangsbedingungen (Zustaende + Nullsensitivitaeten)
XP0    = zeros(nx, np);          % 7 x 8
x0_ext = [x0(:); XP0(:)];        % 7 + 56 = 63 Elemente

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~, X_ext] = ode15s(@(tt,xx) Modell3_woEtOH_XP(tt,xx,u,p,DOTstern), t, x0_ext, opts);

FM    = zeros(np, np);
x_out = zeros(nx, numel(t));
y_out = zeros(ny, numel(t));

for k = 1:numel(t)
    xk         = X_ext(k, 1:nx).';
    x_out(:,k) = xk;
    y_out(:,k) = Modell3_mgl(xk);

    XP_k     = reshape(X_ext(k, nx+1:end), nx, np);
    dhdx     = Modell3_dmgldx(xk);
    dhdtheta = dhdx * XP_k;              % dhdp = 0
    FM       = FM + dhdtheta' * invC * dhdtheta;
end
end
