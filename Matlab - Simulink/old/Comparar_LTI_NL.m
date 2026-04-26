clc

NL_saved = out.NL_saved;
LTI = out.LTI;

% --- Definir intervalo común ---
t0 = max(NL_saved.Time(1), LTI.Time(1));
tf = min(NL_saved.Time(end), LTI.Time(end));

idx = (NL_saved.Time >= t0) & (NL_saved.Time <= tf);
t = NL_saved.Time(idx);

% --- Señales ---
theta_NL = NL_saved.Data(idx,1);

theta_LTI_i = interp1(LTI.Time, LTI.Data(:,1), t, 'linear');

% Eliminar posibles NaN de interpolación
valid = ~isnan(theta_LTI_i);
t = t(valid);
theta_NL = theta_NL(valid);
theta_LTI_i = theta_LTI_i(valid);

% --- Error ---
e_theta = theta_LTI_i - theta_NL;

% =============================
%           FIGURA
% =============================
figure

% --------- (arriba) Error ----------
subplot(2,1,1)
plot(t, e_theta, 'k', 'LineWidth', 1.2)
grid on
ylabel('e_{\theta_m} [rad]')
title('(a) Comparación LTI vs NL - \theta_m')
legend('e(t)','Location','best')
xlim([t(1) t(end)])

% --------- (abajo) Curvas ----------
subplot(2,1,2)
plot(t, theta_NL, 'b', 'LineWidth', 1.2); hold on
plot(t, theta_LTI_i, 'r--', 'LineWidth', 1.2)
grid on
xlabel('Tiempo (s)')
ylabel('\theta_m [rad]')
title('(b) Respuesta \theta_m')
legend('NL','LTI','Location','best')
xlim([t(1) t(end)])