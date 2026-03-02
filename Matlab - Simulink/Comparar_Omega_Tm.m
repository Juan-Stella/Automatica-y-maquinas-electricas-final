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
tiempos_marca = [0 0.3 0.5 0.7 1.2 1.4 1.7 1.9];

for k = 1:length(tiempos_marca)
    [~,i] = min(abs(t - tiempos_marca(k)));

    % Punto grande negro
    plot(omega_NL(i), Tm_NL(i), 'ko', ...
        'MarkerFaceColor','k','MarkerSize',8)

    % Texto del tiempo desplazado
    dx = 30;     % desplazamiento horizontal
    dy = 0.07;   % desplazamiento vertical

    text(omega_NL(i)+dx, Tm_NL(i)+dy, ...
        [num2str(tiempos_marca(k)) ' s'], ...
        'FontSize',10,'FontWeight','bold')
end