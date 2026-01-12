%% VARIABLES DEL SISTEMA
% Proyecto Global Integrador - AyME
% Subsistema mecánico equivalente
% UNCuyo - Ing. Mecatrónica
% Juan Stella - Juan Francisco Huertas

%% -------------------------
%  Carga mecánica
% -------------------------
g = 9.80665;                 % [m/s^2] aceleración de la gravedad

b_l = 0.1;                   % [N*m*s/rad] fricción viscosa en articulación
m   = 1.0;                   % [kg] masa del brazo manipulador
l_cm = 0.25;                 % [m] centro de masa del brazo
J_cm = 0.0208;               % [kg*m^2] inercia del brazo en CM

l_l = 0.50;                  % [m] longitud total del brazo
m_l = 1.5;                   % [kg] masa máxima de carga

% Momento de inercia de la carga
J_l = (m*l_cm^2 + J_cm) + m_l*l_l^2;     % [kg*m^2]
% Coeficiente geométrico del torque gravitacional
k_l = m*l_cm + m_l*l_l;                   % [kg*m]

%% -------------------------
%  Tren de Transmisión
% -------------------------
r = 120;                     % Relación de reducción
Tq_nom = 17;                 % [N*m] Torque nominal salida (regimen continuo o rms)
Tq_max = 45;                 % [N*m] Torque pico salida (corta duración, aceleración)
n_l_nom = 60;                % [rpm] Velocidad nominal salida
omega_l_nom = 6.283;         % [rad/s] Velocidad nominal salida

%% -------------------------
%  Parámetros
% -------------------------
J_m = 1.4e-5;                % [kg*m^2] Inercia del rotor (valor nominal)
b_m = 1.5e-5;                % [N*m*s/rad] Fricción viscosa del motor
Pp = 3;                      % Pares de polos
lambda_m = 0.016;            % [Wb] Flujo magnético equivalente de imanes concatenado por espiras del bobinado de estator   
L_q = 5.8e-3;                % [H] Inductancia de estator (eje en cuadratura)
L_d = 6.6e-3;                % [H] Inductancia de estator (eje directo)
L_ls = 0.8e-3;               % [H] Inductancia de dispersión de estator
R_s = 1.02;                  % [Ω] Resistencia de estator, por fase
alpha_Cu = 3.9e-3;           % [1/°C] Coef. aumento de 𝑅𝑠 con 𝑇𝑠°
C_ts = 0.818;                % [W/°C/s] Capacitancia térmica de estator
R_ts = 146.7;                % [°C/W] Resistencia térmica estator-ambiente
tau_ts_amb = R_ts*C_ts;      % [s] Constante de tiempo térmica

%% -------------------------
%  Especificaciones de operación
% -------------------------
n_m_nom = 6600;              % [rpm] Velocidad nominal del rotor
omega_m_nom = 691.15;        % [rad/s] Velocidad nominal del rotor
V_sl_nom = 30;               % [Vrms] Tensión nominal de línea
V_s_fase_nom = V_sl_nom/sqrt(3); % [Vrms] Tensión nominal de fase
I_s_nom = 0.4;               % [A] Corriente nominal
I_s_max = 2.0;               % [A] Corriente máxima
T_s_max = 115;


%% -------------------------
%  Parámetros equivalentes referidos al eje del motor
% -------------------------
J_eq = J_m + J_l/r^2;        % [kg*m^2] inercia equivalente
b_eq = b_m + b_l/r^2;        % [N*m*s/rad] fricción equivalente






