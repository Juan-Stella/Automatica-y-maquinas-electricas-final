% Variables esperadas:
% outv_abc, outi_abc, outT_m, outomega
clc
r = 120;
[t_v, v_abc]   = extraer_senal(out.v_abc);
[t_i, i_abc]   = extraer_senal(out.i_abc);
[t_T, T_m]     = extraer_senal(out.T_m);
[t_w, omega_m] = extraer_senal(out.omega);

if size(v_abc,2) ~= 3, v_abc = v_abc.'; end
if size(i_abc,2) ~= 3, i_abc = i_abc.'; end

torque_reg  = abs(T_m(end));
torque_pico = max(abs(T_m));

omega_motor_reg  = abs(omega_m(end));
omega_motor_pico = max(abs(omega_m));

omega_caja_reg  = omega_motor_reg/r;
omega_caja_pico = omega_motor_pico/r;

corriente_reg  = max(abs(i_abc(end,:)));
corriente_pico = max(max(abs(i_abc)));

tension_reg  = max(abs(v_abc(end,:)));
tension_pico = max(max(abs(v_abc)));

fprintf('\nRESULTADOS:\n');
fprintf('Torque motor: reg = %.4f, pico = %.4f\n', torque_reg, torque_pico);
fprintf('Velocidad caja: reg = %.4f, pico = %.4f\n', omega_caja_reg, omega_caja_pico);
fprintf('Velocidad motor: reg = %.4f, pico = %.4f\n', omega_motor_reg, omega_motor_pico);
fprintf('Corriente estator: reg = %.4f, pico = %.4f\n', corriente_reg, corriente_pico);
fprintf('Tension de fase: reg = %.4f, pico = %.4f\n', tension_reg, tension_pico);

function [t,y] = extraer_senal(s)
    if isstruct(s)
        t = s.time;
        y = squeeze(s.signals.values);
    elseif isa(s,'timeseries')
        t = s.Time;
        y = squeeze(s.Data);
    else
        error('Formato no soportado');
    end
end