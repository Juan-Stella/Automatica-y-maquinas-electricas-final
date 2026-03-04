%% GRÁFICOS FIELD FORCING/WEAKENING
clear; clc; close all;

load('caso1_nominal.mat');      
load('caso2_forcing.mat');      
load('caso3_weakening.mat');    

%% FIGURA 1: VELOCIDAD
figure;
hold on; grid on;

plot(t_caso1, omega_m_caso1, 'k', 'LineWidth', 2, 'DisplayName', 'v_{ds}^r = 0 V');
plot(t_caso2, omega_m_caso2, 'r--', 'LineWidth', 2, 'DisplayName', 'v_{ds}^r = +1.96 V');
plot(t_caso3, omega_m_caso3, 'b-.', 'LineWidth', 2, 'DisplayName', 'v_{ds}^r = -1.96 V');

xline(0.5, 'k--', 'LineWidth', 1.5);
xlabel('Tiempo [s]', 'FontSize', 12);
ylabel('\omega_m [rad/s]', 'FontSize', 12);
title('Velocidad angular: efecto de field forcing/weakening', 'FontSize', 13);
legend('Location', 'best');

print('velocidad_field_forcing', '-dpng', '-r300');

%% FIGURA 2: TORQUE
figure;
hold on; grid on;

plot(t_caso1, Tm_caso1, 'k', 'LineWidth', 2, 'DisplayName', 'v_{ds}^r = 0 V');
plot(t_caso2, Tm_caso2, 'r--', 'LineWidth', 2, 'DisplayName', 'v_{ds}^r = +1.96 V');
plot(t_caso3, Tm_caso3, 'b-.', 'LineWidth', 2, 'DisplayName', 'v_{ds}^r = -1.96 V');

xline(0.5, 'k--', 'LineWidth', 1.5);
xlabel('Tiempo [s]', 'FontSize', 12);
ylabel('T_m [Nm]', 'FontSize', 12);
title('Torque electromagnético: efecto de field forcing/weakening', 'FontSize', 13);
legend('Location', 'best');

print('torque_field_forcing', '-dpng', '-r300');

disp('Listo.');