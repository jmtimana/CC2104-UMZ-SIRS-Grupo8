%% Grupo 8 - UMZ-SIRS - Runge-Kutta-Fehlberg 4(5) adaptativo (ENTREGA FINAL)
% CC2104 Metodos Numericos. Ejecutar Run con la carpeta matlab/ como directorio.
% RKF45: integrador EXPLICITO adaptativo; no es un solver implicito.
% Salidas:  ../figuras/fig05_paso_adaptativo_rkf45.png
%           ../resultados/MN_Grupo8_resultados_rkf45_tol.csv
% Nota: guardar este archivo en UTF-8 (los textos de las figuras llevan tildes).
clear; clc; close all; format long g

carpetaFiguras='../figuras'; carpetaRes='../resultados';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Parametros y escenarios (los mismos de los archivos 01 y 02)
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=1e-2;
T=100;
condiciones={[0.8;0.2;0;0.8;0.2;0], ...
             [0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original, Z(0)=0','Alternativa, Z(0)=0.1'};
% Tolerancia mixta: evita dividir por valores cercanos a cero.
rtol0=1e-8; atol0=1e-10;      % tolerancia base (tabla de precision)
h0=0.1; hmax=2;               % paso inicial y tope del paso
opt.hmin=1e-10; opt.maxIntentos=200000; opt.tolPos=1e-10; opt.tolCons=1e-8;

%% 2. Ejecucion con la tolerancia base y graficas
for caso=1:numel(condiciones)
    x0=condiciones{caso};
    [t,X,info]=rkf45(x0,p,T,h0,hmax,rtol0,atol0,opt);
    errCons=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
    minP=min(X,[],'all'); nNeg=nnz(X < -opt.tolPos);
    fprintf('\n===== RKF45: %s =====\n',nombres{caso});
    fprintf('Estado=%s | Aceptados=%d, rechazados=%d, eval.f=%d\n', ...
        info.estado,info.aceptados,info.rechazados,info.evalf);
    fprintf('h_min_aceptado=%.5g, h_max_aceptado=%.5g\n',info.hMin,info.hMax);
    fprintf('Error max conservacion=%.3e, minimo=%.3e, negativos=%d\n', ...
        errCons,minP,nNeg);
    fprintf('Estado final [U M Z S I R]:\n'); disp(X(:,end).');

    figure('Color','w','Name',['RKF45 - ' nombres{caso}]);
    subplot(2,2,1);
    plot(t,X(1:3,:)','LineWidth',1.25); grid on;
    xlabel('Tiempo rápido t'); ylabel('Proporción'); title('Capa informativa (RKF45)');
    legend('U','M','Z','Location','best');
    subplot(2,2,2);
    plot(p.eps*t,X(4:6,:)','LineWidth',1.25); grid on;
    xlabel('Tiempo lento \tau=\epsilon t'); ylabel('Proporción'); title('Capa epidemiológica (RKF45)');
    legend('S','I','R','Location','best');
    subplot(2,2,3);
    hh=diff(t); stairs(t(1:end-1),hh,'LineWidth',1.25); grid on;
    xlabel('Tiempo rápido t'); ylabel('Paso aceptado h'); title('Adaptación del paso');
    subplot(2,2,4);
    plot(p.eps*t,X(5,:),'LineWidth',1.4); grid on;
    xlabel('Tiempo lento \tau=\epsilon t'); ylabel('I'); title('Fracción infectada');

    % Figura 05 (solo escenario alternativo): paso aceptado, rechazos y tope hmax.
    if caso==2
        hh=diff(t);
        fig5=figure('Color','w','Name','Figura 05 - Paso adaptativo RKF45');
        % Se excluye el ultimo paso, que se trunca para caer exactamente en T.
        g1=stairs(t(1:end-2),hh(1:end-1),'LineWidth',1.5); hold on;
        leyenda={'Paso aceptado'}; gr=g1;
        if ~isempty(info.tRech)
            g2=plot(info.tRech,info.hRech,'rx','MarkerSize',8,'LineWidth',1.4);
            gr=[gr g2]; leyenda{end+1}='Paso rechazado (h intentado)'; %#ok<AGROW>
        end
        g3=yline(hmax,'--k','LineWidth',1.1); gr=[gr g3]; leyenda{end+1}='Tope h_{max}';
        grid on; xlabel('Tiempo rápido t'); ylabel('Paso h');
        title(sprintf('Adaptación del paso en RKF45 (\\epsilon=%g, tol. relativa %g, Z(0)=0,1)',p.eps,rtol0));
        legend(gr,leyenda,'Location','southeast');
        % fig05 ya no se usa en el paper (se eliminó el Anexo C); no se exporta.
        % exportgraphics(fig5,fullfile(carpetaFiguras,'fig05_paso_adaptativo_rkf45.png'),'Resolution',300);
    end
end

%% 3. Barrido de tolerancias (estudio de convergencia correcto para RKF45)
% En un metodo adaptativo se varia la TOLERANCIA, no el paso inicial h0.
% Referencia: ode45 con tolerancias estrictas (ver tambien el archivo 09).
optRef=odeset('RelTol',1e-12,'AbsTol',1e-13,'MaxStep',0.01);
rtols=10.^(-4:-1:-10);
filas={};
for caso=1:numel(condiciones)
    x0=condiciones{caso};
    solRef=ode45(@(tt,x)modelo(x,p),[0 T],x0,optRef);
    for k=1:numel(rtols)
        rt=rtols(k); at=1e-2*rt;
        reloj=tic;
        [t,X,info]=rkf45(x0,p,T,h0,hmax,rt,at,opt);
        seg=toc(reloj);
        if strcmp(info.estado,'OK')
            err=max(abs(X-deval(solRef,min(t,T))),[],'all');
            cons=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
            minimo=min(X,[],'all');
        else
            err=NaN; cons=NaN; minimo=NaN;
        end
        filas(end+1,:)={caso,rt,at,err,cons,minimo,info.aceptados,info.rechazados, ...
            info.evalf,info.hMin,info.hMax,seg,info.estado}; %#ok<AGROW>
        fprintf('Caso %d | RelTol=%.0e | error=%.3e | aceptados=%d | rechazados=%d | evalf=%d\n', ...
            caso,rt,err,info.aceptados,info.rechazados,info.evalf);
    end
end
tabla=cell2table(filas,'VariableNames',{'Caso','RelTol','AbsTol','ErrorMaxRef', ...
    'ErrorConservacion','Minimo','Aceptados','Rechazados','EvalF','hMinAceptado', ...
    'hMaxAceptado','Segundos','Estado'});
writetable(tabla,fullfile(carpetaRes,'MN_Grupo8_resultados_rkf45_tol.csv'));
fprintf('\nGuardado: MN_Grupo8_resultados_rkf45_tol.csv\n');
fprintf('RKF45 COMPLETADO. Los errores son controles numericos, no errores exactos.\n');

%% FUNCIONES LOCALES
function [t,X,info]=rkf45(x0,p,T,h0,hmax,rtol,atol,opt)
    cap=1024; tb=zeros(1,cap); Xb=zeros(6,cap); n=1; Xb(:,1)=x0;
    h=h0; intentos=0;
    info=struct('estado','OK','aceptados',0,'rechazados',0,'evalf',0,'rechPos',0, ...
        'hMin',inf,'hMax',0,'tRech',[],'hRech',[]);
    while T-tb(n)>1e-12*max(1,T)
        intentos=intentos+1;
        if intentos>opt.maxIntentos || h<opt.hmin
            info.estado='FALLO (paso minimo o demasiados intentos)'; break
        end
        h=min([h,hmax,T-tb(n)]);
        x=Xb(:,n);
        k1=modelo(x,p);
        k2=modelo(x+h*(1/4*k1),p);
        k3=modelo(x+h*(3/32*k1+9/32*k2),p);
        k4=modelo(x+h*(1932/2197*k1-7200/2197*k2+7296/2197*k3),p);
        k5=modelo(x+h*(439/216*k1-8*k2+3680/513*k3-845/4104*k4),p);
        k6=modelo(x+h*(-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5),p);
        info.evalf=info.evalf+6;
        x4=x+h*(25/216*k1+1408/2565*k3+2197/4104*k4-1/5*k5);
        x5=x+h*(16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6);
        esc=atol+rtol.*max(abs(x),abs(x5));
        er=max(abs(x5-x4)./esc);
        if ~isfinite(er), er=inf; end
        valido=all(isfinite(x5)) && all(x5>=-opt.tolPos) && 1-p.K*x5(5)>0 && 1+x5(3)>0;
        if er<=1 && valido
            n=n+1; if n>cap, cap=2*cap; tb(cap)=0; Xb(6,cap)=0; end
            tb(n)=tb(n-1)+h; Xb(:,n)=x5;   % solucion de orden 5
            info.aceptados=info.aceptados+1;
            info.hMin=min(info.hMin,h); info.hMax=max(info.hMax,h);
        else
            info.rechazados=info.rechazados+1;
            if ~valido, info.rechPos=info.rechPos+1; end
            info.tRech(end+1)=tb(n); info.hRech(end+1)=h; %#ok<AGROW>
        end
        if er==0, fac=2; else, fac=min(2,max(0.2,0.9*(1/er)^(1/5))); end
        if ~valido, fac=min(fac,0.5); end
        h=h*fac;
    end
    t=tb(1:n); X=Xb(:,1:n);
end

function dx=modelo(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5); R=x(6);
    if 1-p.K*I<=0 || 1+Z<=0, dx=NaN(6,1); return, end
    a=p.b2/(1-p.K*I); b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end
