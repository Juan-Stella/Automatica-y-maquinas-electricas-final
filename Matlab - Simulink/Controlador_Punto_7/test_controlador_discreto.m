%% SCRIPT DE PRUEBA: test_controlador_discreto.m
% Este script permite probar (correr) la función controlador_discreto sin
% errores de "Not enough input arguments", ya que le pasaremos parámetros.

% 1. Limpiamos variables, consola y persistentes previos (importante)
clear all; close all; clc;

% 2. Ejecutamos los scripts base del proyecto para cargar las variables
% Asegurarnos de que las variables como r, J_eq, b_a, etc., estén en memoria.
% (Descomenta la siguiente línea si necesitas que cargue directo todo el workspace)
% run('../variables.m'); 
% run('../Polos_con_Jeq.m'); 

% - Valores manuales de prueba por defecto si no están en workspace:
r = 120;
Ts = 1e-3; % Asumiendo 1 milisegundo de muestreo
K = 1;
J_eq = 1.4e-5 + (0.0208 + 1.5*0.5^2)/(120^2); 
b_a = 0.05;
K_sa = 10;
K_sia = 50;
k_theta2 = 0.8;
k_omega2 = 1.2;
k_i = 100;

% 3. Configuración de Entradas Ficticias (para una prueba rápida)
q_ref_test     = 10.0; % Consigna de posición de carga [rad]
q_dot_ref_test = 0.0;  % Consigna de velocidad de carga [rad/s]
theta_m_test   = 0.0;  % Posición leída del motor [rad] (ej: encoder)
sel_pos_test   = 1.0;  % 1: Habilitar control de Posición

% 4. Llamado a la función
% Primera iteración k=1
try
    [T_m_ref(1), theta_m_obs(1), omega_m_obs(1)] = controlador_discreto(...
        q_ref_test, q_dot_ref_test, theta_m_test, sel_pos_test, ...
        Ts, K, r, J_eq, b_a, K_sa, K_sia, k_theta2, k_omega2, k_i);
    
    disp('La iteración 1 funcionó correctamente. T_ref:');
    disp(T_m_ref(1));
    
    % Segunda iteración k=2
    [T_m_ref(2), theta_m_obs(2), omega_m_obs(2)] = controlador_discreto(...
        q_ref_test, q_dot_ref_test, theta_m_test, sel_pos_test, ...
        Ts, K, r, J_eq, b_a, K_sa, K_sia, k_theta2, k_omega2, k_i);
        
    disp('La iteración 2 funcionó correctamente. La memoria funciona.');
catch ME
    disp('Hubo un error evaluando:');
    disp(ME.message);
end
