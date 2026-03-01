clc
omega_NL_saved;
omega_LTI_saved;
% --- Definir intervalo común ---
t0 = max(omega_NL_saved.Time(1), omega_LTI_saved.Time(1));
tf = min(omega_NL_saved.Time(end), omega_LTI_saved.Time(end));

idx = (omega_NL_saved.Time >= t0) & (omega_NL_saved.Time <= tf);
t = omega_NL_saved.Time(idx);

% --- Señales ---
theta_NL = omega_NL_saved.Data(idx,1);

theta_omega_LTI_saved_i = interp1(omega_LTI_saved.Time, omega_LTI_saved.Data(:,1), t, 'linear');

% Eliminar posibles NaN de interpolación
valid = ~isnan(theta_omega_LTI_saved_i);
t = t(valid);
theta_NL = theta_NL(valid);
theta_omega_LTI_saved_i = theta_omega_LTI_saved_i(valid);

% --- Error ---
e_theta = theta_omega_LTI_saved_i - theta_NL;

% =============================
%           FIGURA
% =============================
figure

% --------- (arriba) Error ----------
subplot(2,1,1)
plot(t, e_theta, 'k', 'LineWidth', 1.2)
grid on
ylabel('e_{\theta_m} [rad]')
title('(a) Comparación omega_LTI_saved vs NL - \theta_m')
legend('e(t)','Location','best')
xlim([t(1) t(end)])

% --------- (abajo) Curvas ----------
subplot(2,1,2)
plot(t, theta_NL, 'b', 'LineWidth', 1.2); hold on
plot(t, theta_omega_LTI_saved_i, 'r--', 'LineWidth', 1.2)
grid on
xlabel('Tiempo (s)')
ylabel('\theta_m [rad]')
title('(b) Respuesta \theta_m')
legend('NL','omega_LTI_saved','Location','best')
xlim([t(1) t(end)])