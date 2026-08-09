% check_jacobians.m
%
% Prueft die symbolisch erzeugten Ableitungen gegen zentrale finite
% Differenzen des RHS. Muss nach jedem Lauf von gen_jacobians.m ausgefuehrt
% werden -- eine falsche Jacobi-Matrix faellt sonst nirgends auf, sie
% verfaelscht nur still die FIM.
%
% Erwartung: Abweichungen in der Groessenordnung 1e-9 relativ zu den
% Eintraegen. Die beiden Normen am Ende pruefen gezielt die Spalten, die in
% einer aelteren Version identisch null waren (KPh bzw. mPh) -- damals war
% der Phosphat-Term in gen_jacobians.m schlicht vergessen worden.

clear; clc;
projectRoot = pwd;
addpath(fullfile(projectRoot,'..','Modelle'),'-begin');
addpath(fullfile(projectRoot,'..','utils'),'-begin');

load(fullfile(projectRoot,'Daten','Daten_Processed', ...
              'Processed_FedBatch_Modell3_MultiExp.mat'));

u_t      = TrainSet(1).u;
DOTstern = max(TrainSet(1).O2.y);
t_t      = 25;                        % Zeitpunkt mit aktivem Feed

% Auswertepunkt: realistischer Zustand, keine Nullen (sonst sind viele
% Ableitungen trivial null und der Test sagt nichts aus).
xt = [10; 50; 20; 5; 10; 0.2; 80];    % V, mX, mGlc, mAm, mPh, mB, DOT
p = [0.3430 4.29 0.1435 0.0852 0.0837 0.0241 1.0 246.0375 0.01 0.01];

%% df/dx ------------------------------------------------------------------
A_num = zeros(7,7);
for j = 1:7
    h  = 1e-6 * max(abs(xt(j)),1);
    xp = xt;  xp(j) = xp(j)+h;
    xm = xt;  xm(j) = xm(j)-h;
    A_num(:,j) = (Modell3_woEtOH_10p(t_t,xp,u_t,p,DOTstern) - ...
                  Modell3_woEtOH_10p(t_t,xm,u_t,p,DOTstern)) / (2*h);
end
A_sym = Modell3_dfdx_10p(xt, p, DOTstern);
fprintf('max |dfdx_sym - dfdx_num| = %.3e\n', max(abs(A_sym(:)-A_num(:))));

%% df/dp ------------------------------------------------------------------
B_num = zeros(7,numel(p));
for j = 1:numel(p)
    h  = 1e-6 * max(abs(p(j)),1);
    pp = p;  pp(j) = pp(j)+h;
    pm = p;  pm(j) = pm(j)-h;
    B_num(:,j) = (Modell3_woEtOH_10p(t_t,xt,u_t,pp,DOTstern) - ...
                  Modell3_woEtOH_10p(t_t,xt,u_t,pm,DOTstern)) / (2*h);
end
B_sym = Modell3_dfdp_10p(xt, p, DOTstern);
fprintf('max |dfdp_sym - dfdp_num| = %.3e\n', max(abs(B_sym(:)-B_num(:))));

%% Gezielte Kontrolle: Phosphat-Term vorhanden? --------------------------
fprintf('||dfdp(:,10)|| = %.3e  (KPh, muss > 0 sein)\n', norm(B_sym(:,10)));
fprintf('||dfdx(:,5)||  = %.3e  (mPh, muss > 0 sein)\n', norm(A_sym(:,5)));

%% Messgleichung ----------------------------------------------------------
% dmgldx wandelt Zustands- in Ausgangssensitivitaeten um (Massen ->
% Konzentrationen) und steckt damit in derselben Kette wie oben.
J_num = zeros(6,7);
for j = 1:7
    h  = 1e-6 * max(abs(xt(j)),1);
    xp = xt;  xp(j) = xp(j)+h;
    xm = xt;  xm(j) = xm(j)-h;
    J_num(:,j) = (Modell3_mgl(xp) - Modell3_mgl(xm)) / (2*h);
end
fprintf('max |dmgldx_sym - dmgldx_num| = %.3e\n', ...
        max(abs(Modell3_dmgldx(xt) - J_num), [], 'all'));