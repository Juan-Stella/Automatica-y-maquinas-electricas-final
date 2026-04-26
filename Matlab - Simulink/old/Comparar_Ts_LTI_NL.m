clc
Ts_NL_saved = out.Ts_NL_saved;
Ts_LTI_saved = out.Ts_LTI_saved;
% --- Definir intervalo común ---
t0 = max(Ts_NL_saved.Time(1), Ts_LTI_saved.Time(1));
tf = min(Ts_NL_saved.Time(end), Ts_LTI_saved.Time(end));

idx = (Ts_NL_saved.Time >= t0) & (Ts_NL_saved.Time <= tf);
t = Ts_NL_saved.Time(idx);

% --- Señales ---
theta_NL = Ts_NL_saved.Data(idx,1);

theta_Ts_LTI_saved_i = interp1(Ts_LTI_saved.Time, Ts_LTI_saved.Data(:,1), t, 'linear');

% Eliminar posibles NaN de interpolación
valid = ~isnan(theta_Ts_LTI_saved_i);
t = t(valid);
theta_NL = theta_NL(valid);
theta_Ts_LTI_saved_i = theta_Ts_LTI_saved_i(valid);

% --- Error ---
e_theta = theta_Ts_LTI_saved_i - theta_NL;

% =============================
%           FIGURA
% =============================
figure

% --------- (arriba) Error ----------
subplot(2,1,1)
plot(t, e_theta, 'k', 'LineWidth', 1.2)
grid on
ylabel('e_{T_s} [°C]')
title('(a) Comparación T_s LTI vs NL')
legend('e(t)','Location','best')
xlim([t(1) t(end)])

% --------- (abajo) Curvas ----------
subplot(2,1,2)
plot(t, theta_NL, 'b', 'LineWidth', 1.2); hold on
plot(t, theta_Ts_LTI_saved_i, 'r--', 'LineWidth', 1.2)
grid on
xlabel('Tiempo (s)')
ylabel('T_s [°C]')
title('(b) Respuesta T_s')
legend('T_s NL','T_s LTI','Location','best')
xlim([t(1) t(end)])