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

T_s_fix = 40;                 % [°C] Temperatura fija para este análisis
T_ref_Rs = 20;                % [°C] Referencia de resistencia (Aleo)
R_s_fix = R_s_ref*(1 + alpha_Cu*(T_s_fix - T_ref_Rs));

% --- Valores nominales (ya definidos en tu script)
J_eq_nom = J_eq;
b_eq_nom = b_eq;

% --- Valores mínimos (ej: m_l = 0 kg, b_l = 0)
m_l_min = 0;                  % [kg]
b_l_min = 0;                  % [N*m*s/rad]

J_l_min = (m*l_cm^2 + J_cm) + m_l_min*l_l^2;
J_eq_min = J_m + J_l_min/r^2;

b_eq_min = b_m + b_l_min/r^2;

% --- Valores máximos (ej: m_l = m_l, b_l = b_l)
m_l_max = m_l;                % [kg] (tu valor máximo definido)
b_l_max = b_l;                % [N*m*s/rad] (tu valor máximo definido)

J_l_max = (m*l_cm^2 + J_cm) + m_l_max*l_l^2;
J_eq_max = J_m + J_l_max/r^2;

b_eq_max = b_m + b_l_max/r^2;

% Guardamos en vectores para iterar
J_eq_vec = [J_eq_nom, J_eq_min, J_eq_max];
b_eq_vec = [b_eq_nom, b_eq_min, b_eq_max];

% Resultados
p1_v = zeros(1,3);
p2_v = zeros(1,3);
p3_v = zeros(1,3);
cero_v = zeros(1,3);
wn_v = zeros(1,3);
zeta_v = zeros(1,3);

for k = 1:3
    J_eq_k = J_eq_vec(k);
    b_eq_k = b_eq_vec(k);

    A = L_q*b_eq_k + J_eq_k*R_s_fix;
    B = J_eq_k*L_q;
    C = R_s_fix*b_eq_k + (3/2)*Pp^2*lambda_m^2;

    Delta = A^2 - 4*B*C;

    p2_v(k) = (-A + sqrt(Delta))/(2*B);
    p3_v(k) = (-A - sqrt(Delta))/(2*B);

    cero_v(k) = -R_s_fix/L_q;

    wn_v(k) = abs(p2_v(k));
    zeta_v(k) = -real(p2_v(k))/wn_v(k);
end

% Etiquetas filas como en el informe
Condicion = ["Nominal"; "Minimos"; "Maximos"];

Tabla_param = table( ...
    Condicion, ...
    J_eq_vec.', ...
    b_eq_vec.', ...
    p1_v.', ...
    p2_v.', ...
    p3_v.', ...
    cero_v.', ...
    wn_v.', ...
    zeta_v.', ...
    'VariableNames', ...
    {'Condicion','J_eq','b_eq','P1','P2','P3','Cero','wn_rad_s','zeta'} );

disp(Tabla_param)
writetable(Tabla_param,'Tabla_Polos_Jeq_beq.xlsx')
