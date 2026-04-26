clc

omega_NL_saved;
Tm_NL_saved;

t = omega_NL_saved.Time;
omega_NL = omega_NL_saved.Data(:,1);
Tm_NL    = Tm_NL_saved.Data(:,1);

% Recorte 0–2 s
idx = (t>=0 & t<=2);
t = t(idx);
omega_NL = omega_NL(idx);
Tm_NL = Tm_NL(idx);

figure
plot(omega_NL, Tm_NL,'b','LineWidth',1.3); hold on
grid on
xlabel('\omega_m [rad/s]')
ylabel('T_m [N·m]')
title('Curva torque vs velocidad (NL)')
xline(0,'k'); yline(0,'k');

% --------- TIEMPOS QUE QUERÉS MARCAR ----------
tiempos_marca = [0 0.1001 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9];

% Arrays para guardar valores
omega_vals = zeros(length(tiempos_marca),1);
Tm_vals    = zeros(length(tiempos_marca),1);

for k = 1:length(tiempos_marca)

    [~,i] = min(abs(t - tiempos_marca(k)));

    omega_vals(k) = omega_NL(i);
    Tm_vals(k)    = Tm_NL(i);

    % Punto negro
    plot(omega_NL(i), Tm_NL(i), 'ko', ...
        'MarkerFaceColor','k','MarkerSize',8)

    % Texto desplazado
    dx = 30;
    dy = 0.07;

    text(omega_NL(i)+dx, Tm_NL(i)+dy, ...
        [num2str(tiempos_marca(k)) ' s'], ...
        'FontSize',10,'FontWeight','bold')
end

% ===============================
%        MOSTRAR RESULTADOS
% ===============================

tabla_resultados = table(tiempos_marca.', omega_vals, Tm_vals, ...
    'VariableNames', {'Tiempo_s','Omega_rad_s','Tm_Nm'});

disp(tabla_resultados)

% ===============================
%     OPCIONAL: Exportar Excel
% ===============================

writetable(tabla_resultados,'valores_curva_NL.xlsx');