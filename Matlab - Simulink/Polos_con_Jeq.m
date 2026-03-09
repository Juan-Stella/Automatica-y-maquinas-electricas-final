%% Polos en el plano complejo - formato similar al usado antes

zeta = 0.75;
omega_pos = 800;      % [rad/s]
polo_corriente = -5000;

%% -------------------------
%  Polos a lazo abierto
% -------------------------

T_s = T_s_ref;
R_s = R_s_ref*(1 + alpha_Cu*(T_s - T_s_ref));

A = L_q*b_eq + J_eq*R_s;
B = J_eq*L_q;
C = R_s*b_eq + (3/2)*Pp^2*lambda_m^2;
Delta = A^2 - 4*B*C;

p1_ab = 0;
p2_ab = (-A + sqrt(Delta)) / (2*B);
p3_ab = (-A - sqrt(Delta)) / (2*B);

polos_ab = [p1_ab; p2_ab; p3_ab];

%% -------------------------
%  Ganancias PID
% -------------------------

b_a   = J_eq * (2*zeta + 1) * omega_pos
K_sa  = J_eq * (2*zeta + 1) * omega_pos^2
K_sia = J_eq * omega_pos^3

%% -------------------------
%  Polos PID - nominal
% -------------------------

P_nom = [J_eq b_a K_sa K_sia];
polos_pid_nom = roots(P_nom);

%% -------------------------
%  Polos PID - máximo
% -------------------------

P_max = [J_eq_max b_a K_sa K_sia];
polos_pid_max = roots(P_max);

%% -------------------------
%  Gráfico
% -------------------------

figure
hold on
grid on
box on

plot(real(polos_ab), imag(polos_ab), 'xk', 'LineWidth', 1.5, 'MarkerSize', 7)
plot(real(polos_pid_nom), imag(polos_pid_nom), 'xr', 'LineWidth', 1.5, 'MarkerSize', 7)
plot(real(polos_pid_max), imag(polos_pid_max), 'xb', 'LineWidth', 1.5, 'MarkerSize', 7)
plot(real(polo_corriente), 0, 'xg', 'LineWidth', 1.5, 'MarkerSize', 7)

xlabel('Parte Real')
ylabel('Parte Imaginaria')
title('Ubicación de Polos en el Plano Complejo')

legend('Polos Lazo Abierto', ...
       'Polos PID - valores nominales/minimos', ...
       'Polos PID - valores máximos', ...
       'Polo lazo corriente', ...
       'Location', 'northeast')

xlim([-5200 100])
ylim([-1600 1600])