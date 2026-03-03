clc
omega_NL_saved; iqs_NL_saved;

t  = omega_NL_saved.Time;
w  = omega_NL_saved.Data(:,1);
iq = iqs_NL_saved.Data(:,1);

idx = (t>=0 & t<=2);
t=t(idx); w=w(idx); iq=iq(idx);

eventos = [0.1 0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.9];  % ordenados

settle_pct = 0.01;  % ±1%
dt_pre = 0.01; dt_fin = 0.02;

% umbral solo para overshoot (%)
Dymin_w  = 1.0;     % rad/s
Dymin_iq = 0.05;    % A

tabla_omega = tablaIdx(t,w, eventoss(eventos),settle_pct,dt_pre,dt_fin,Dymin_w);
tabla_iqs   = tablaIdx(t,iq,eventoss(eventos),settle_pct,dt_pre,dt_fin,Dymin_iq);

tabla_omega.Properties.VariableNames = {'Perturbacion_s','RiseFall_ms','Settling_ms','OverUnder_pct','ValorEstab'};
tabla_iqs.Properties.VariableNames   = {'Perturbacion_s','RiseFall_ms','Settling_ms','OverUnder_pct','ValorEstab'};

writetable(tabla_omega,'Resultados_transitorios.xlsx','Sheet','Velocidad');
writetable(tabla_iqs,'Resultados_transitorios.xlsx','Sheet','Corriente');

% ===== FUNCIONES =====
function e = eventoss(e), e = sort(e); end

function T = tablaIdx(t,y,eventos,settle_pct,dt_pre,dt_fin,Dymin)
n=length(eventos);
Rise_ms=zeros(n,1); Settle_ms=zeros(n,1); OS_pct=zeros(n,1); Yf=zeros(n,1);

for k=1:n
    t0=eventos(k);
    if k<n, t1=eventos(k+1); else, t1=t(end); end

    idxw = (t>=t0 & t<=t1);
    tt = t(idxw); yy = y(idxw);

    y0 = mean(y(t>=t0-dt_pre & t<t0));
    yf = mean(y(t>=t1-dt_fin & t<=t1)); Yf(k)=yf;
    Dy = yf - y0;

    % ---- Rise/Fall 10-90 (si no alcanza, usar ventana completa) ----
    y10 = y0 + 0.1*Dy;  y90 = y0 + 0.9*Dy;
    if Dy >= 0
        i10=find(yy>=y10,1,'first'); i90=find(yy>=y90,1,'first');
    else
        i10=find(yy<=y90,1,'first'); i90=find(yy<=y10,1,'first');
    end
    if isempty(i10) || isempty(i90)
        Rise_ms(k) = (t1-t0)*1000;   % no llegó dentro de la ventana
    else
        Rise_ms(k) = abs(tt(i90)-tt(i10))*1000;
    end

    % ---- Settling ±1% (si no asienta, usar ventana completa) ----
    band = settle_pct*max(abs(Dy),Dymin);   % banda mínima
    err = abs(yy - yf);
    Settle_ms(k) = (t1-t0)*1000;            % por defecto: no asentó
    for i=1:length(tt)
        if all(err(i:end)<=band)
            Settle_ms(k) = (tt(i)-t0)*1000;
            break
        end
    end

    % ---- Over/Under (%) (si Dy chico, 0%) ----
    if abs(Dy) < Dymin
        OS_pct(k) = 0;
    else
        if Dy > 0, OS_pct(k)=max(0,(max(yy)-yf)/abs(Dy))*100;
        else,      OS_pct(k)=max(0,(yf-min(yy))/abs(Dy))*100; end
    end
end

T = table(eventos(:),Rise_ms,Settle_ms,OS_pct,Yf);
end