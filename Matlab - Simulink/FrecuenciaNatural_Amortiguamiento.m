clc

% ===== DATOS =====
omega_NL_saved; iqs_NL_saved;
t = omega_NL_saved.Time;
omega = omega_NL_saved.Data(:,1);
iqs   = iqs_NL_saved.Data(:,1);

% Recorte 0–2 s
idx = (t>=0 & t<=2);
t = t(idx); omega = omega(idx); iqs = iqs(idx);

% ===== EVENTOS =====
eventos = [0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9];

% ===== PARAMETROS =====
settle_pct = 0.01; rise_low=0.1; rise_high=0.9; dt_pre=0.01; dt_fin=0.02;

% ===== TABLAS (se crean ACÁ) =====
tabla_omega = calcTabla(t, omega, eventos, settle_pct, rise_low, rise_high, dt_pre, dt_fin);
tabla_iqs   = calcTabla(t, iqs,   eventos, settle_pct, rise_low, rise_high, dt_pre, dt_fin);

% Renombrar columnas como Aleo
tabla_omega.Properties.VariableNames = {'Perturbacion_s','RiseTimeFallTime_ms','SettlingTime_ms','OverUnder_pct','ValorEstablecimiento'};
tabla_iqs.Properties.VariableNames   = {'Perturbacion_s','RiseTimeFallTime_ms','SettlingTime_ms','OverUnder_pct','ValorEstablecimiento'};

% Mostrar
disp('=== Resultados para velocidad (omega_m) ==='); disp(tabla_omega)
disp('=== Resultados para corriente (i_qs) ===');    disp(tabla_iqs)

% Export a Excel (tu estilo)
writetable(tabla_omega,'Resultados_velocidad.xlsx');
writetable(tabla_iqs,'Resultados_corriente.xlsx');

% =======================
%        FUNCION
% =======================
function T = calcTabla(t,y,eventos,settle_pct,rise_low,rise_high,dt_pre,dt_fin)

n = length(eventos);
Rise_ms=nan(n,1); Settle_ms=nan(n,1); OS_pct=nan(n,1); Yf=nan(n,1);

for k=1:n
    t0 = eventos(k);
    if k<n, t1 = eventos(k+1); else, t1 = t(end); end

    idxw=(t>=t0)&(t<=t1); tt=t(idxw); yy=y(idxw);

    y0 = mean(y((t>=(t0-dt_pre))&(t<t0)));
    yf = mean(y((t>=(t1-dt_fin))&(t<=t1))); Yf(k)=yf;

    Dy = yf-y0; if abs(Dy)<1e-9, continue, end

    yn=(yy-y0)/Dy;
    i10=find(yn>=rise_low,1,'first'); i90=find(yn>=rise_high,1,'first');
    if ~isempty(i10)&&~isempty(i90), Rise_ms(k)=(tt(i90)-tt(i10))*1000; end

    band=settle_pct*abs(Dy); err=abs(yy-yf);
    for i=1:length(tt)
        if all(err(i:end)<=band), Settle_ms(k)=(tt(i)-t0)*1000; break, end
    end

    if Dy>0, OS_pct(k)=max(0,(max(yy)-yf)/abs(Dy))*100;
    else,    OS_pct(k)=max(0,(yf-min(yy))/abs(Dy))*100; end
end

T = table(eventos(:), Rise_ms, Settle_ms, OS_pct, Yf);
end