clc
Tm_NL_saved;
Tm_LTI_saved;
% --- Definir intervalo común ---
t0 = max(Tm_NL_saved.Time(1), Tm_LTI_saved.Time(1));
tf = min(Tm_NL_saved.Time(end), Tm_LTI_saved.Time(end));

idx = (Tm_NL_saved.Time >= t0) & (Tm_NL_saved.Time <= tf);
t = Tm_NL_saved.Time(idx);

% --- Señales ---
theta_NL = Tm_NL_saved.Data(idx,1);

theta_Tm_LTI_saved_i = interp1(Tm_LTI_saved.Time, Tm_LTI_saved.Data(:,1), t, 'linear');

% Eliminar posibles NaN de interpolación
valid = ~isnan(theta_Tm_LTI_saved_i);
t = t(valid);
theta_NL = theta_NL(valid);
theta_Tm_LTI_saved_i = theta_Tm_LTI_saved_i(valid);

% --- Error ---
e_theta = theta_Tm_LTI_saved_i - theta_NL;

% =============================
%           FIGURA
% =============================
figure

% --------- (arriba) Error ----------
subplot(2,1,1)
plot(t, e_theta, 'k', 'LineWidth', 1.2)
grid on
ylabel('e_{T_m} [N*m]')
title('(a) Comparación T_m LTI vs NL')
legend('e(t)','Location','best')
xlim([t(1) t(end)])

% --------- (abajo) Curvas ----------
subplot(2,1,2)
plot(t, theta_NL, 'b', 'LineWidth', 1.2); hold on
plot(t, theta_Tm_LTI_saved_i, 'r--', 'LineWidth', 1.2)
grid on
xlabel('Tiempo (s)')
ylabel('T_m [N*m]')
title('(b) Respuesta T_m')
legend('NL','T_m LTI','Location','best')
xlim([t(1) t(end)])