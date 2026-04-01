function [T_m_ref, theta_m_obs, omega_m_obs] = controlador_discreto(q_ref, q_dot_ref, theta_m, sel_pos, Ts, K, r, J_eq, b_a, K_sa, K_sia, k_theta2, k_omega2, k_i)
% CONTROLADOR_DISCRETO  Controlador programado en Texto Estructurado (Punto 7).
% Reemplaza el modelo en bloques del controlador discreto completo 
% compuesto por: Set-Point (Fig 92), PID externo (Fig 93) y Observador Aumentado (Fig 94).
%
% Entradas:
%   q_ref     : Referencia de posición de carga [rad]
%   q_dot_ref : Referencia de velocidad de carga [rad/s]
%   theta_m   : Posición medida del motor [rad] (retroalimentación)
%   sel_pos   : Selector de modo (1: Regula posición, 0: Regula velocidad)
%   [... 10 Parámetros ...] : Variables de diseño provenientes del workspace.
%
% Salidas:
%   T_m_ref     : Torque de control de referencia para el motor [Nm]
%   theta_m_obs : Posición del motor observada [rad]
%   omega_m_obs : Velocidad del motor observada [rad/s]
%
% Documentado con una estructura típica de desarrollo de Microcontrolador/DSP.
% Integrantes: Juan Stella - Juan Francisco Huertas

    % =====================================================================
    % 1. DECLARACIÓN DE ESTADOS (VARIABLES PERSISTENTES)
    % =====================================================================
    % Simulando el entorno de un microcontrolador (variables static o globales),
    % mantienen la "memoria" del sistema de lazo discreto entre ejecuciones (k-1).
    
    persistent is_initialized
    
    % Constantes de precálculo (ahorran ciclos de ALU en el bucle principal)
    persistent C_r C_invTs C_Tustin C_invJeq
    persistent C_ba C_Ksa C_Ksia
    persistent C_ktheta2 C_komega2 C_ki 
    
    % a) Estados del generador de Set-Point (Figura 92)
    persistent sp_err_prev
    persistent sp_tust_u_prev sp_tust_y_prev
    
    % b) Estados del lazo PID Externo (Figura 93)
    persistent pid_tust1_u_prev pid_tust1_y_prev
    persistent pid_tust2_u_prev pid_tust2_y_prev
    
    % c) Estados del Observador Aumentado (Figura 94)
    persistent obs_z1 obs_z2 obs_z3
    persistent obs_tust1_u_prev obs_tust1_y_prev
    persistent obs_tust2_u_prev obs_tust2_y_prev
    persistent obs_tust3_u_prev obs_tust3_y_prev
    
    % =====================================================================
    % 2. RUTINA DE INICIALIZACIÓN Y PRECÁLCULO (Se ejecuta una sola vez)
    % =====================================================================
    if isempty(is_initialized)
        is_initialized = true;
        
        % -> Precálculo de operaciones fijas:
        C_r       = r;
        C_invTs   = 1 / Ts;
        C_Tustin  = K * Ts / 2;
        C_invJeq  = 1 / J_eq;
        
        C_ba      = b_a;
        C_Ksa     = K_sa;
        C_Ksia    = K_sia;
        
        C_ktheta2 = k_theta2;
        C_komega2 = k_omega2;
        C_ki      = k_i;
        
        % -> Set a cero de registros tipo retardo (condiciones iniciales nulas)
        sp_err_prev      = 0;
        sp_tust_u_prev   = 0; sp_tust_y_prev   = 0;
        
        pid_tust1_u_prev = 0; pid_tust1_y_prev = 0;
        pid_tust2_u_prev = 0; pid_tust2_y_prev = 0;
        
        obs_z1           = 0; obs_z2           = 0; obs_z3           = 0;
        obs_tust1_u_prev = 0; obs_tust1_y_prev = 0;
        obs_tust2_u_prev = 0; obs_tust2_y_prev = 0;
        obs_tust3_u_prev = 0; obs_tust3_y_prev = 0;
    end
    
    % =====================================================================
    % 3. CÁLCULO DE SALIDAS AL INSTANTE ACTUAL "k" (Bucle principal)
    % =====================================================================
    
    % --- A. EXTRAER SALIDAS DEL OBSERVADOR AUMENTADO ---
    % Como existen registros de retardo puro (1/z) previos a los integradores, 
    % los lazos algebraicos están desacoplados. Las salidas se evalúan directo.
    
    % Rama 1: Integrador de Error del Observador
    [int_err, obs_tust1_u_prev, obs_tust1_y_prev] = func_Tustin(obs_z1, obs_tust1_u_prev, obs_tust1_y_prev, C_Tustin);
    
    % Rama 2: Obtención de Velocidad (omega_m_obs)
    [omega_m_obs, obs_tust2_u_prev, obs_tust2_y_prev] = func_Tustin(obs_z2, obs_tust2_u_prev, obs_tust2_y_prev, C_Tustin);
    
    % Rama 3: Obtención de Posición (theta_m_obs)
    [theta_m_obs, obs_tust3_u_prev, obs_tust3_y_prev] = func_Tustin(obs_z3, obs_tust3_u_prev, obs_tust3_y_prev, C_Tustin);
    
    
    % --- B. CÁLCULO DEL SET-POINT ---
    % Adaptación para referir las señales desde la carga al eje del motor
    q_ref_motor     = C_r * q_ref;
    q_dot_ref_motor = C_r * q_dot_ref;
    
    % Rama Derivadora: de posición a velocidad
    q_ref_deriv = (q_ref_motor - sp_err_prev) * C_invTs;
    sp_err_prev = q_ref_motor;  % Se actualiza retardo z^-1
    
    % Rama Integradora: de velocidad a posición (no explotada por este PID, pero presente en diagrama)
    [~, sp_tust_u_prev, sp_tust_y_prev] = func_Tustin(q_dot_ref_motor, sp_tust_u_prev, sp_tust_y_prev, C_Tustin);
    
    % Selector o llave de consignas (1: Control por SetPoint de Posicion)
    if sel_pos > 0.5
        omega_m_ref = q_ref_deriv;
    else
        omega_m_ref = q_dot_ref_motor;
    end
    
    
    % --- C. CÁLCULO DEL PID EXTERNO ---
    err_omega = omega_m_ref - omega_m_obs;
    
    % Primer Integrador (Integra error de velocidad -> Error de Posición)
    [pid_x1, pid_tust1_u_prev, pid_tust1_y_prev] = func_Tustin(err_omega, pid_tust1_u_prev, pid_tust1_y_prev, C_Tustin);
    
    % Segundo Integrador (Integra Error de Posición -> Error de Posición Acumulado)
    [pid_x2, pid_tust2_u_prev, pid_tust2_y_prev] = func_Tustin(pid_x1, pid_tust2_u_prev, pid_tust2_y_prev, C_Tustin);
    
    % Torque de Referencia Demandado al Motor (Salida)
    T_m_ref = (C_ba * err_omega) + (C_Ksa * pid_x1) + (C_Ksia * pid_x2);
    
    
    % --- D. RETROALIMENTACIÓN DE ESTADOS DEL OBSERVADOR ---
    % Se evalúan los nodos pre-retardos que se usarán en "k+1"
    err_obs = theta_m - theta_m_obs;
    k_i_out = C_ki * int_err;
    
    % Nodos sumadores
    accel_node = (T_m_ref * C_invJeq) + (C_komega2 * err_obs) + k_i_out; % Aceleración
    vel_node   = omega_m_obs + (C_ktheta2 * err_obs);                    % Velocidad corregida
    
    % Actualización de registros de retardo puros (bloques 1/z)
    obs_z1 = err_obs;
    obs_z2 = accel_node;
    obs_z3 = vel_node;
    
end

% =========================================================================
% FUNCIONES ESCALARES DESAGREGADAS
% =========================================================================

function [y, u_prev, y_prev] = func_Tustin(u, u_prev, y_prev, C_Tustin)
% FUNC_TUSTIN: Integrador discreto mediante método Trapezoide (Tustin).
% Abstrae la operación para no repetir código, facilitando la lectura.
% 
% Ecuación de recurrencia en el tiempo:
% y[k] = y[k-1] + (K * Ts / 2) * (u[k] + u[k-1])
%
% Entradas:
%   u        : Valor de entrada en 'k'
%   u_prev   : Valor de entrada en 'k-1'
%   y_prev   : Salida acumulada en 'k-1'
%   C_Tustin : Ganancia precalculada integradora
%
% Salidas:
%   y        : Nueva salida en 'k'
%   u_prev   : Valor para registrar en delay
%   y_prev   : Valor para registrar en delay

    y = y_prev + C_Tustin * (u + u_prev);
    
    % Desplazamiento para el próximo instante (Overwriting)
    u_prev = u;
    y_prev = y;
end
