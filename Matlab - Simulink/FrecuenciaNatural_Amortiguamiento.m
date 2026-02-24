clc

%% ---------------------------------------------------------
%  CUADRO: Comparación de wn y zeta bajo diferentes condiciones
% ---------------------------------------------------------

% === Resistencias a 40°C y 115°C ===
T_ref_Rs = 20;

R_s_40  = R_s_ref*(1 + alpha_Cu*(40  - T_ref_Rs));
R_s_115 = R_s_ref*(1 + alpha_Cu*(115 - T_ref_Rs));

% === J_eq y b_eq nominales ===
J_eq_nom = J_eq;
b_eq_nom = b_eq;

% === J_eq y b_eq mínimos ===
m_l_min = 0;                  % [kg]
b_l_min = 0.07;               % [N*m*s/rad]

J_l_min  = (m*l_cm^2 + J_cm) + m_l_min*l_l^2;
J_eq_min = J_m + J_l_min/r^2;
b_eq_min = b_m + b_l_min/r^2;

% === J_eq y b_eq máximos ===
m_l_max = 1.5;                % [kg]
b_l_max = 0.13;               % [N*m*s/rad]

J_l_max  = (m*l_cm^2 + J_cm) + m_l_max*l_l^2;
J_eq_max = J_m + J_l_max/r^2;
b_eq_max = b_m + b_l_max/r^2;

% === Armar condiciones (6 filas) ===
Condicion = strings(6,1);
wn_val    = zeros(6,1);
zeta_val  = zeros(6,1);

% ---- 1) Manteniendo valores equivalentes nominales ----
i = 1;
R_s = R_s_40;  J_eq = J_eq_nom;  b_eq = b_eq_nom;
Condicion(i) = "R_s = " + string(R_s) + " (40°C)";
wn_val(i)   = sqrt( (R_s*b_eq + (3/2)*Pp^2*(lambda_m^2)) / (J_eq*L_q) );
zeta_val(i) = (L_q*b_eq + R_s*J_eq) / (2*J_eq*L_q*wn_val(i));

i = 2;
R_s = R_s_115; J_eq = J_eq_nom;  b_eq = b_eq_nom;
Condicion(i) = "R_s = " + string(R_s) + " (115°C)";
wn_val(i)   = sqrt( (R_s*b_eq + (3/2)*Pp^2*(lambda_m^2)) / (J_eq*L_q) );
zeta_val(i) = (L_q*b_eq + R_s*J_eq) / (2*J_eq*L_q*wn_val(i));

% ---- 2) Manteniendo R_s a 40°C y modificando J_eq y b_eq ----
i = 3;
R_s = R_s_40;  J_eq = J_eq_min;  b_eq = b_eq_min;
Condicion(i) = "J_{eq}, b_{eq} mínimos (R_s 40°C)";
wn_val(i)   = sqrt( (R_s*b_eq + (3/2)*Pp^2*(lambda_m^2)) / (J_eq*L_q) );
zeta_val(i) = (L_q*b_eq + R_s*J_eq) / (2*J_eq*L_q*wn_val(i));

i = 4;
R_s = R_s_40;  J_eq = J_eq_max;  b_eq = b_eq_max;
Condicion(i) = "J_{eq}, b_{eq} máximos (R_s 40°C)";
wn_val(i)   = sqrt( (R_s*b_eq + (3/2)*Pp^2*(lambda_m^2)) / (J_eq*L_q) );
zeta_val(i) = (L_q*b_eq + R_s*J_eq) / (2*J_eq*L_q*wn_val(i));

% ---- 3) Manteniendo R_s a 115°C y modificando J_eq y b_eq ----
i = 5;
R_s = R_s_115; J_eq = J_eq_min;  b_eq = b_eq_min;
Condicion(i) = "J_{eq}, b_{eq} mínimos (R_s 115°C)";
wn_val(i)   = sqrt( (R_s*b_eq + (3/2)*Pp^2*(lambda_m^2)) / (J_eq*L_q) );
zeta_val(i) = (L_q*b_eq + R_s*J_eq) / (2*J_eq*L_q*wn_val(i));

i = 6;
R_s = R_s_115; J_eq = J_eq_max;  b_eq = b_eq_max;
Condicion(i) = "J_{eq}, b_{eq} máximos (R_s 115°C)";
wn_val(i)   = sqrt( (R_s*b_eq + (3/2)*Pp^2*(lambda_m^2)) / (J_eq*L_q) );
zeta_val(i) = (L_q*b_eq + R_s*J_eq) / (2*J_eq*L_q*wn_val(i));

% === Tabla final ===
Cuadro3 = table(Condicion, wn_val, zeta_val, 'VariableNames', {'Condicion','wn_rad_s','zeta'});
disp(Cuadro3);

% Export a Excel
writetable(Cuadro3,'Cuadro3_wn_zeta.xlsx');

% Impresión con 4 decimales
fprintf('\n%-45s | %12s | %10s\n','Condición','wn [rad/s]','zeta');
fprintf('%s\n', repmat('-',1,75));
for i = 1:height(Cuadro3)
    fprintf('%-45s | %12.4f | %10.4f\n', Cuadro3.Condicion(i), Cuadro3.wn_rad_s(i), Cuadro3.zeta(i));
end