clc

% ===== Señales (timeseries) =====
omega_NL_saved;
iqs_NL_saved;

t  = omega_NL_saved.Time;
w  = omega_NL_saved.Data(:,1);
iq = iqs_NL_saved.Data(:,1);

% ===== Recorte 0–2 s =====
idx = (t>=0 & t<=2);
t=t(idx); w=w(idx); iq=iq(idx);

% ===== Eventos (como Aleo en tabla) =====
eventos = [0.1 0.3 0.5 0.7 0.9 1.1 1.5 1.7 1.9];

% ===== Parámetros =====
settle_pct = 0.01;   % ±1%
dt0 = 0.01;
dtf = 0.02;

% ===== Tablas =====
tabla_w  = calcTablaAleo(t,w,  eventos, settle_pct, dt0, dtf, 1.0);
tabla_iq = calcTablaAleo(t,iq, eventos, settle_pct, dt0, dtf, 0.05);

tabla_w.Properties.VariableNames  = {'Perturbacion_s','RiseTimeFallTime_ms','SettlingTime_ms','OvershootUndershoot_pct','ValorEstablecimiento_rad_s'};
tabla_iq.Properties.VariableNames = {'Perturbacion_s','RiseTimeFallTime_ms','SettlingTime_ms','OvershootUndershoot_pct','ValorEstablecimiento_A'};

disp('=== Resultados para velocidad (omega_m) ==='); disp(tabla_w)
disp('=== Resultados para corriente (i_qs) ===');    disp(tabla_iq)

% ===== Export =====
archivo = 'Resultados_Aleo.xlsx';
writetable(tabla_w,  archivo, 'Sheet','Velocidad');
writetable(tabla_iq, archivo, 'Sheet','Corriente');

% ==========================================================
function T = calcTablaAleo(t,y,eventos,settle_pct,dt0,dtf,Dymin)

n = numel(eventos);
Rise_ms   = nan(n,1);
Settle_ms = nan(n,1);
OS_pct    = nan(n,1);
Yf        = nan(n,1);

for k=1:n
    t0 = eventos(k);
    if k<n, t1 = eventos(k+1); else, t1 = t(end); end

    idxw = (t>=t0) & (t<=t1);
    tt = t(idxw);
    yy = y(idxw);
    if numel(tt)<3, continue, end

    idx0 = (t >= (t0-dt0)) & (t < t0);
    idxf = (t > (t1-dtf)) & (t <= t1);
    if ~any(idx0) || ~any(idxf), continue, end

    y0 = mean(y(idx0));
    yf = mean(y(idxf));
    Dy = yf - y0;

    % Rise 10–90
    y10 = y0 + 0.1*Dy;
    y90 = y0 + 0.9*Dy;

    if Dy >= 0
        i10 = find(yy >= y10, 1, 'first');
        i90 = find(yy >= y90, 1, 'first');
    else
        i10 = find(yy <= y90, 1, 'first');
        i90 = find(yy <= y10, 1, 'first');
    end

    if isempty(i10) || isempty(i90) || i90==i10
        Rise_ms(k) = (t1 - t0)*1000;
    else
        Rise_ms(k) = abs(tt(i90) - tt(i10))*1000;
    end

    % Settling ±1%
    band = settle_pct * max(abs(Dy), Dymin);
    err = abs(yy - yf);

    Settle_ms(k) = (t1 - t0)*1000;
    for i=1:numel(tt)
        if all(err(i:end) <= band)
            Settle_ms(k) = (tt(i) - t0)*1000;
            break
        end
    end

    % Overshoot/Undershoot %
    if abs(Dy) < Dymin
        OS_pct(k) = 0;
    else
        if Dy > 0
            OS_pct(k) = max(0,(max(yy)-yf)/abs(Dy))*100;
        else
            OS_pct(k) = max(0,(yf-min(yy))/abs(Dy))*100;
        end
    end

    Yf(k) = yf;
end

% ===== REDONDEO Y FORMATO FIJO 4 DECIMALES =====
Rise_ms   = string(compose("%.4f", Rise_ms));
Settle_ms = string(compose("%.4f", Settle_ms));
OS_pct    = string(compose("%.4f", OS_pct));
Yf        = string(compose("%.4f", Yf));

eventos_str = string(compose("%.4f", eventos(:)));

T = table(eventos_str, Rise_ms, Settle_ms, OS_pct, Yf);
end