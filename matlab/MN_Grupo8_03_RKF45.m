%% Grupo 8 - UMZ-SIRS - Runge-Kutta-Fehlberg 4(5) adaptativo
% CC2104 Metodos Numericos. Ejecutar Run en MATLAB.
% RKF45: integrador EXPLICITO adaptativo; no es un solver implicito.
clear; clc; close all; format long g

%% Exportacion automatica para el paper
carpetaFiguras = 'figuras';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end

%% 1. Parametros y escenarios ya usados en el archivo 02
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=1e-2;
T=100;
condiciones={[0.8;0.2;0;0.8;0.2;0], ...
             [0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original, Z(0)=0','Alternativa, Z(0)=0.1'};
% Tolerancia mixta: evita division por valores cercanos a cero.
opt.atol=1e-10; opt.rtol=1e-8;
opt.h0=0.1; opt.hmin=1e-9; opt.hmax=2;
opt.maxPasos=200000; opt.tolCons=1e-8; opt.tolPos=1e-10;

%% 2. Ejecutar RKF45 con control de error local
for caso=1:numel(condiciones)
    x0=condiciones{caso};
    [t,X,info]=integrar_rkf45(x0,p,T,opt);

    errU=max(abs(sum(X(1:3,:),1)-1));
    errS=max(abs(sum(X(4:6,:),1)-1));
    errCons=max(errU,errS);
    minP=min(X,[],'all');
    nNeg=nnz(X < -opt.tolPos);
    fprintf('\n===== RKF45: %s =====\n',nombres{caso});
    fprintf('Aceptados=%d, rechazados=%d, eval.f=%d\n', ...
        info.aceptados,info.rechazados,info.evalf);
    fprintf('h_min_aceptado=%.5g, h_max_aceptado=%.5g\n', ...
        info.hminAceptado,info.hmaxAceptado);
    fprintf('Error max conservacion=%.3e, minimo=%.3e, negativos=%d\n', ...
        errCons,minP,nNeg);
    fprintf('Conservacion: %s | Positividad: %s\n', ...
        estado(errCons < opt.tolCons),estado(nNeg==0));
    fprintf('Estado final [U M Z S I R]:\n'); disp(X(:,end).');

    figure('Color','w','Name',['RKF45 - ' nombres{caso}]);
    subplot(2,2,1);
    plot(t,X(1:3,:)','LineWidth',1.25); grid on;
    xlabel('t rapido'); ylabel('Proporcion'); title('Informacion: RKF45');
    legend('U','M','Z','Location','best');

    subplot(2,2,2);
    plot(p.eps*t,X(4:6,:)','LineWidth',1.25); grid on;
    xlabel('tau = epsilon*t'); ylabel('Proporcion'); title('Epidemia: RKF45');
    legend('S','I','R','Location','best');

    subplot(2,2,3);
    stairs(t(1:end-1),diff(t),'LineWidth',1.25); grid on;
    xlabel('t rapido'); ylabel('h aceptado'); title('Adaptacion del paso');

    subplot(2,2,4);
    plot(p.eps*t,X(5,:),'LineWidth',1.4); grid on;
    xlabel('tau = epsilon*t'); ylabel('I'); title('Fraccion infectada');

    % Figura 05: solo el escenario alternativo, en una figura clara para el paper.
    if caso==2
        figPaso=figure('Color','w','Name','Figura 05 - Paso adaptativo RKF45');
        stairs(t(1:end-1),diff(t),'LineWidth',1.5);
        grid on; xlabel('Tiempo rapido t'); ylabel('Paso aceptado h');
        title('Adaptacion del paso RKF45: Z(0)=0.1');
        exportgraphics(figPaso,fullfile(carpetaFiguras, ...
            'fig05_paso_adaptativo_rkf45.png'),'Resolution',300);
        fprintf('Figura exportada: fig05_paso_adaptativo_rkf45.png\n');
    end
end
fprintf('\nRKF45 COMPLETADO. Los errores son controles numericos, no errores exactos.\n');

%% FUNCIONES LOCALES
function [t,X,info]=integrar_rkf45(x0,p,T,opt)
    t=zeros(1,opt.maxPasos+1);
    X=zeros(6,opt.maxPasos+1); X(:,1)=x0;
    h=opt.h0; n=1; rechazados=0; evalf=0;
    hMin=inf; hMax=0; intentos=0;
    while t(n)<T
        intentos=intentos+1;
        if intentos>opt.maxPasos
            error('RKF45 supero el limite de intentos.');
        end
        h=min([h,opt.hmax,T-t(n)]);
        if h<opt.hmin && T-t(n)>opt.hmin
            error('RKF45: h inferior a hmin.');
        end
        [x4,x5,nf]=paso_rkf45(X(:,n),h,p);
        evalf=evalf+nf;
        escala=opt.atol+opt.rtol.*max(abs(X(:,n)),abs(x5));
        errorEsc=max(abs(x5-x4)./escala);
        % El control de dominio rechaza pasos fisicamente inadmisibles.
        dominioValido=all(isfinite(x5)) && all(x5>=-opt.tolPos) ...
            && 1-p.K*x5(5)>0 && 1+x5(3)>0;
        if errorEsc<=1 && dominioValido
            n=n+1;
            t(n)=t(n-1)+h;
            X(:,n)=x5;  % solucion de orden 5
            hMin=min(hMin,h); hMax=max(hMax,h);
        else
            rechazados=rechazados+1;
        end
        if ~isfinite(errorEsc), errorEsc=inf; end
        if errorEsc==0
            fac=2;
        else
            fac=min(2,max(0.2,0.9*(1/errorEsc)^(1/5)));
        end
        if ~dominioValido, fac=min(fac,0.5); end
        h=h*fac;
    end
    t=t(1:n); X=X(:,1:n);
    info.aceptados=n-1; info.rechazados=rechazados;
    info.evalf=evalf; info.hminAceptado=hMin; info.hmaxAceptado=hMax;
end

function [x4,x5,nf]=paso_rkf45(x,h,p)
    % Coeficientes clasicos de Fehlberg (ordenes 4 y 5).
    k1=f(x,p);
    k2=f(x+h*(1/4*k1),p);
    k3=f(x+h*(3/32*k1+9/32*k2),p);
    k4=f(x+h*(1932/2197*k1-7200/2197*k2+7296/2197*k3),p);
    k5=f(x+h*(439/216*k1-8*k2+3680/513*k3-845/4104*k4),p);
    k6=f(x+h*(-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5),p);
    x4=x+h*(25/216*k1+1408/2565*k3+2197/4104*k4-1/5*k5);
    x5=x+h*(16/135*k1+6656/12825*k3+28561/56430*k4-9/50*k5+2/55*k6);
    nf=6;
end

function dx=f(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5); R=x(6);
    if 1-p.K*I<=0 || 1+Z<=0, error('Denominador no positivo.'); end
    a=p.b2/(1-p.K*I); b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end

function s=estado(ok)
    if ok, s='OK'; else, s='NO CUMPLE'; end
end
