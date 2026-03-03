%% CÁLCULO Y GRÁFICO DEL ÁNGULO DE DESFASAJE δ(t)


t = out.v_a.Time;                   
v_a = out.v_a.Data;          
v_b = out.v_b.Data;         
v_c = out.v_c.Data;        
theta_m = out.Theta_NL.Data;    
Pp = 3;  

% Fórmulas de Clarke
v_alpha = (2/3) * (v_a - 0.5*v_b - 0.5*v_c);
v_beta = (2/3) * (sqrt(3)/2 * v_b - sqrt(3)/2 * v_c);

%calculo sin corrección de saltos
theta_ev_raw = atan2(v_beta, v_alpha);

%corregir saltos
theta_ev = unwrap(theta_ev_raw);

theta_r = Pp * theta_m;

%desfasaje:
delta = theta_ev - theta_r;

% Grafico
figure('Position', [100, 100, 900, 500]);
plot(t, rad2deg(delta), 'Color', [0.6 0 0.6], 'LineWidth', 2);
grid on;
xlabel('t [s]', 'FontSize', 12);
ylabel('\delta(t) [°]', 'FontSize', 12);
title('Evolución temporal del ángulo de desfasaje electromagnético', ...
      'FontSize', 13, 'FontWeight', 'bold');

% Línea de referencia en δ=0
yline(0, 'k--', 'LineWidth', 1.5);

% Configuración de ejes
set(gca, 'FontSize', 11);
xlim([0 max(t)]);

% Guardar figura
print(gcf, 'angulo_desfasaje', '-dpng', '-r300');

% est
fprintf('\n=== ESTADÍSTICAS DEL ÁNGULO δ(t) ===\n');
fprintf('δ máximo:  %+.2f° (%+.4f rad)\n', max(rad2deg(delta)), max(delta));
fprintf('δ mínimo:  %+.2f° (%+.4f rad)\n', min(rad2deg(delta)), min(delta));
fprintf('δ promedio: %+.2f° (%+.4f rad)\n', mean(rad2deg(delta)), mean(delta));