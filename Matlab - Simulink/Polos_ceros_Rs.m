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