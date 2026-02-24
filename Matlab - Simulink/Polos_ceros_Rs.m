clc

T_s_ref = 20;
T_s_interval = 40:10:115;
n = numel(T_s_interval);

p2 = zeros(1,n);
p3 = zeros(1,n);
ceros = zeros(1,n);
R_s = zeros(1,n);
for i = 1:n
    T_s = T_s_interval(i);

    R_s(i) = R_s_ref*(1 + alpha_Cu*(T_s - T_s_ref));

    A = L_q*b_eq + J_eq*R_s(i);
    B = J_eq*L_q;
    C = R_s(i)*b_eq + (3/2)*Pp^2*lambda_m^2;

    Delta = A^2 - 4*B*C;

    p2(i) = (-A + sqrt(Delta)) / (2*B);
    p3(i) = (-A - sqrt(Delta)) / (2*B);

    ceros(i) = -R_s(i)/L_q;


 
end

wn = abs(p2);                % omega_n
zeta = -real(p2)./wn;         % amortiguamiento

figure
plot(real(p2), imag(p2), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5); hold on
plot(real(p3), imag(p3), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(real(ceros), imag(ceros), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5)

% Flecha polos s2
quiver(real(p2(1)), imag(p2(1)), ...
       real(p2(end)) - real(p2(1)), ...
       imag(p2(end)) - imag(p2(1)), ...
       0, 'r', 'LineWidth', 1.5)

% Flecha polos s3
quiver(real(p3(1)), imag(p3(1)), ...
       real(p3(end)) - real(p3(1)), ...
       imag(p3(end)) - imag(p3(1)), ...
       0, 'b', 'LineWidth', 1.5)

% Flecha ceros
quiver(real(ceros(1)), imag(ceros(1)), ...
       real(ceros(end)) - real(ceros(1)), ...
       imag(ceros(end)) - imag(ceros(1)), ...
       0, 'k', 'LineWidth', 1.5)

% Etiquetas para polos s2
dy = 0; % ajuste vertical
text(real(p2(1))+2, imag(p2(1))+dy,  ' 40°C')
text(real(p2(end))+2,imag(p2(end))+dy,' 115°C')

% Etiquetas para polos s3
dy = 5; % ajuste vertical
text(real(p3(1)), imag(p3(1))+dy,  ' 40°C')
text(real(p3(end)),imag(p3(end))+dy,' 115°C')

% Etiquetas para ceros
dy = 7; % ajuste vertical
text(real(ceros(1))-4,  imag(ceros(1))+dy,  ' 40°C')
text(real(ceros(end))-4,imag(ceros(end))+dy,' 115°C')

grid on
xlabel('Parte Real')
ylabel('Parte Imaginaria')
title('Variación de Polos y Ceros con T_s')
legend('Polo s_2','Polo s_3','Cero','Location','best')

axis equal


% Polo s1 (siempre 0)
p1 = zeros(1,n);

% Crear tabla
Tabla = table( ...
    R_s.', ...
    p1.', ...
    p2.', ...
    p3.', ...
    ceros.', ...
    wn.', ...
    zeta.', ...
    'VariableNames', ...
    {'Rs_Ohm','P1','P2','P3','Ceros','wn_rad_s','zeta'} );

disp(Tabla)
writetable(Tabla,'Tabla_Polos_Temperatura.xlsx')

%% ---------------------------------------------------------
%  Variación de parámetros mecánicos: J_eq y b_eq (T_s = 40°C)
% ---------------------------------------------------------

T_s_fix = 40;                
T_ref_Rs = 20;                
R_s_fix = R_s_ref*(1 + alpha_Cu*(T_s_fix - T_ref_Rs));

% --- Valores nominales (ya definidos en tu script)
J_eq_nom = J_eq;
b_eq_nom = b_eq;

% --- Valores mínimos
m_l_min = 0;                  % [kg]
b_l_min = 0.07;               % [N*m*s/rad]

J_l_min = (m*l_cm^2 + J_cm) + m_l_min*l_l^2;
J_eq_min = J_m + J_l_min/r^2;

b_eq_min = b_m + b_l_min/r^2;

% --- Valores máximo
m_l_max = 1.5;                % [kg] 
b_l_max = 0.13;               % [N*m*s/rad] 

J_l_max = (m*l_cm^2 + J_cm) + m_l_max*l_l^2;
J_eq_max = J_m + J_l_max/r^2;

N = 10; 

J_eq_sweep = linspace(J_eq_min, J_eq_max, N);
b_eq_sweep = linspace(b_eq_min, b_eq_max, N);

p1_s = zeros(1,N);
p2_s = zeros(1,N);
p3_s = zeros(1,N);
cero_s = zeros(1,N);
wn_s = zeros(1,N);
zeta_s = zeros(1,N);

% Resultados
p1_v = zeros(1,N);
p2_v = zeros(1,N);
p3_v = zeros(1,N);
cero_v = zeros(1,N);
wn_v = zeros(1,N);
zeta_v = zeros(1,N);

for k = 1:N
    Jk = J_eq_sweep(k);
    bk = b_eq_sweep(k);

    A = L_q*bk + Jk*R_s_fix;
    B = Jk*L_q;
    C = R_s_fix*bk + (3/2)*Pp^2*lambda_m^2;

    Delta = A^2 - 4*B*C;

    p2_s(k) = (-A + sqrt(Delta))/(2*B);
    p3_s(k) = (-A - sqrt(Delta))/(2*B);

    cero_s(k) = -R_s_fix/L_q;

    wn_s(k) = abs(p2_s(k));
    zeta_s(k) = -real(p2_s(k))/wn_s(k);
end

% Export tabla del barrido (esto sí te da "muchos puntos")
Tabla_param_sweep = table( ...
    J_eq_sweep.', ...
    b_eq_sweep.', ...
    p1_s.', ...
    p2_s.', ...
    p3_s.', ...
    cero_s.', ...
    wn_s.', ...
    zeta_s.', ...
    'VariableNames', {'J_eq','b_eq','P1','P2','P3','Cero','wn_rad_s','zeta'} );

writetable(Tabla_param_sweep,'Tabla_Sweep_Jeq_beq.xlsx');

% ---------------- Figura ----------------
figure
plot(real(p1_s), imag(p1_s), 'mx', 'MarkerSize', 8, 'LineWidth', 1.5); hold on
plot(real(p2_s), imag(p2_s), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(real(p3_s), imag(p3_s), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5)
plot(real(cero_s), imag(cero_s), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5)

% Flechas al costado (sin superponer)
xoff2 = max(real(p2_s)) + 10; % desplazamiento a la derecha (ajustable)
y2a = imag(p2_s(1));
y2b = imag(p2_s(end));
quiver(xoff2, y2a, 0, (y2b - y2a), 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.7)
text(xoff2+2, (y2a+y2b)/2, 'De J_{eq} y b_{eq} mínimo a máximo', ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle')

xoff3 = max(real(p3_s)) + 10;
y3a = imag(p3_s(1));
y3b = imag(p3_s(end));
quiver(xoff3, y3a, 0, (y3b - y3a), 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.7)
text(xoff3+2, (y3a+y3b)/2, 'De J_{eq} y b_{eq} mínimo a máximo', ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle')

grid on
xlabel('Parte Real')
ylabel('Parte Imaginaria')
title('Ubicación de Polos y Ceros en el Plano Complejo')
legend('Polos s1','Polos s2','Polos s3','Cero','Location','best')
axis equal