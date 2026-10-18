%% GRUPO 8 - Comparacion cuantitativa de integradores UMZ-SIRS (CORREGIDO)
% MATLAB R2019b+ | Guardar como Live Script (.mlx) desde Live Editor.
% Solucion de referencia: ode45 con tolerancias estrictas (NO exacta).
% Se miden errores de trayectoria en nodos comunes y tiempo de ejecucion.
clear; clc; close all; format long g

%% Exportacion automatica para el paper
carpetaFiguras = 'figuras';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end

%% 1. Parametros y controles
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=1e-2;
T=100;
hs=[0.4 0.2 0.1 0.05]; % dividen exactamente T
condiciones={[0.8;0.2;0;0.8;0.2;0],[0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original Z0=0','Alternativa Z0=0.1'};
metodos={'Euler explicito','RK4','RKF45','Euler atras'};
optRef=odeset('RelTol',1e-12,'AbsTol',1e-13,'MaxStep',0.01);
optRK.rtol=1e-8; optRK.atol=1e-10; optRK.hmax=2; optRK.hmin=1e-10;
optImp.tolG=1e-10; optImp.tolGS=1e-11;
optImp.maxNewton=50; optImp.maxGS=50;

%% 2. Comparacion en cada condicion inicial
for caso=1:2
    x0=condiciones{caso};
    fprintf('\n================ %s ================\n',nombres{caso});
    tic;
    solRef=ode45(@(t,x)modelo(x,p),[0 T],x0,optRef);
    tref=solRef.x;
    fprintf('Referencia ode45: %.3f s, %d nodos internos\n',toc,numel(tref));
    % Los metodos de paso fijo se comparan en sus nodos; RKF45 en los propios.
    resultados={};
    for ih=1:numel(hs)
        h=hs(ih); N=round(T/h); t=(0:N)*h;
        % deval usa la interpolacion continua interna de ode45.
        Yref=deval(solRef,t);
        for im=1:4
            nombre=metodos{im};
            try
                reloj=tic;
                if im==1
                    [X,feval]=expl_fijo(x0,p,h,N,1); nSteps=N; nRej=0; itN=0; itGS=0;
                elseif im==2
                    [X,feval]=expl_fijo(x0,p,h,N,4); nSteps=N; nRej=0; itN=0; itGS=0;
                elseif im==3
                    [ta,Xa,feval,nRej]=rkf45(x0,p,T,h,optRK);
                    X=Xa; % Para RKF45 evaluar en sus PROPIOS nodos, sin pchip externo.
                    Ycomparar=deval(solRef,ta);
                    nSteps=numel(ta)-1; itN=0; itGS=0;
                else
                    [X,feval,itN,itGS]=impl_fijo(x0,p,h,N,optImp);
                    nSteps=N; nRej=0;
                end
                segundos=toc(reloj);
                if im~=3, Ycomparar=Yref; end
                err=max(abs(X-Ycomparar),[],'all');
                conserv=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
                minimo=min(X,[],'all');
                resultados(end+1,:)={h,nombre,err,conserv,minimo,nSteps,nRej,feval,itN,itGS,segundos,'OK'}; %#ok<AGROW>
            catch ME
                resultados(end+1,:)={h,nombre,NaN,NaN,NaN,0,0,0,0,0,NaN,['FALLO: ' ME.message]}; %#ok<AGROW>
            end
        end
    end
    tabla=cell2table(resultados,'VariableNames', ...
        {'h','Metodo','ErrorMaxRef','ErrorConservacion','Minimo','Pasos','Rechazados','EvalF','IterNewton','IterGS','Segundos','Estado'});
    disp(tabla);
    nombreCSV=sprintf('MN_Grupo8_resultados_caso_%d.csv',caso);
    writetable(tabla,nombreCSV);
    fprintf('Tabla guardada: %s\n',nombreCSV);

    figure('Color','w','Name',['Comparacion - ' nombres{caso}]);
    hold on;
    for im=1:4
        filas=strcmp(tabla.Metodo,metodos{im}) & strcmp(tabla.Estado,'OK');
        if im==3, continue; end % RKF45 adapta h: variar h0 no es estudio de orden.
        if any(filas)
            loglog(tabla.h(filas),max(tabla.ErrorMaxRef(filas),eps),'o-','LineWidth',1.4,'DisplayName',metodos{im});
        end
    end
    grid on; xlabel('Paso inicial o fijo h'); ylabel('Error maximo vs referencia');
    title(['Convergencia (paso fijo): ' nombres{caso}]); legend('Location','best');
    set(gca,'XScale','log','YScale','log');
    % Figura 01: caso alternativo, con escala logaritmica para distinguir RK4.
    if caso==2
        exportgraphics(gcf,fullfile(carpetaFiguras, ...
            'fig01_convergencia_alternativa.png'),'Resolution',300);
        fprintf('Figura exportada: fig01_convergencia_alternativa.png\n');
    end
end
fprintf('\nCOMPARACION COMPLETADA. Revisar posibles FALLOS y discrepancias.\n');
fprintf('El tiempo de ejecucion depende de MATLAB, equipo y carga de trabajo.\n');
fprintf('RKF45 se compara en sus nodos aceptados: no se interpola su trayectoria.\n');
fprintf('Para estudiar convergencia de RKF45 debe variarse RelTol y AbsTol, no h0.\n');

%% Funciones locales (no modificar funciones de plantillas oficiales, si las hubiera)
function dx=modelo(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5); R=x(6);
    a=p.b2/(1-p.K*I); b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end

function [X,nf]=expl_fijo(x,p,h,N,orden)
    X=zeros(6,N+1); X(:,1)=x; nf=0;
    for n=1:N
        if orden==1
            X(:,n+1)=X(:,n)+h*modelo(X(:,n),p); nf=nf+1;
        else
            a=modelo(X(:,n),p);
            b=modelo(X(:,n)+h*a/2,p);
            c=modelo(X(:,n)+h*b/2,p);
            d=modelo(X(:,n)+h*c,p);
            X(:,n+1)=X(:,n)+h*(a+2*b+2*c+d)/6; nf=nf+4;
        end
        if any(~isfinite(X(:,n+1))) || (1-p.K*X(5,n+1)<=0) || (1+X(3,n+1)<=0)
            error('Trayectoria fuera de dominio');
        end
    end
end

function [ts,Xs,nf,rech]=rkf45(x,p,T,h,o)
    ts=0; Xs=x; nf=0; rech=0; intento=0;
    while ts(end)<T-1e-12
        intento=intento+1;
        if intento>200000 || h<o.hmin, error('RKF45: demasiados intentos/paso pequeno'); end
        h=min([h,o.hmax,T-ts(end)]);
        x=Xs(:,end);
        k1=modelo(x,p);
        k2=modelo(x+h*(k1/4),p);
        k3=modelo(x+h*(3*k1/32+9*k2/32),p);
        k4=modelo(x+h*(1932*k1/2197-7200*k2/2197+7296*k3/2197),p);
        k5=modelo(x+h*(439*k1/216-8*k2+3680*k3/513-845*k4/4104),p);
        k6=modelo(x+h*(-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40),p);
        nf=nf+6;
        y4=x+h*(25*k1/216+1408*k3/2565+2197*k4/4104-k5/5);
        y5=x+h*(16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55);
        sc=o.atol+o.rtol.*max(abs(x),abs(y5));
        er=max(abs(y5-y4)./sc);
        valido=all(isfinite(y5)) && all(y5>=-1e-10) && 1-p.K*y5(5)>0 && 1+y5(3)>0;
        if er<=1 && valido
            ts(end+1)=ts(end)+h; %#ok<AGROW>
            Xs(:,end+1)=y5; %#ok<AGROW>
        else
            rech=rech+1;
        end
        if er==0, fac=2; else, fac=min(2,max(0.2,0.9*er^(-1/5))); end
        if ~valido, fac=min(fac,0.5); end
        h=h*fac;
    end
end

function [X,nf,itN,itGS]=impl_fijo(x,p,h,N,o)
    X=zeros(6,N+1); X(:,1)=x; nf=0; itN=0; itGS=0;
    for n=1:N
        xn=X(:,n); y=xn; ok=false;
        for k=1:o.maxNewton
            F=modelo(y,p); nf=nf+1;
            G=y-xn-h*F;
            if norm(G,inf)<o.tolG, ok=true; break; end
            W=eye(6)-h*jac_f(y,p);
            if any(abs(diag(W))<1e-12), error('Diagonal casi nula'); end
            D=diag(diag(W)); L=tril(W,-1); U=triu(W,1);
            if max(abs(eig(-(D+L)\U)))>=1, error('Gauss-Seidel no converge'); end
            delta=zeros(6,1); gsOK=false;
            for j=1:o.maxGS
                ant=delta;
                for i=1:6
                    delta(i)=(-G(i)-W(i,1:i-1)*delta(1:i-1)-W(i,i+1:6)*ant(i+1:6))/W(i,i);
                end
                itGS=itGS+1;
                if norm(W*delta+G,inf)<o.tolGS, gsOK=true; break; end
            end
            if ~gsOK, error('GS no alcanza residuo'); end
            alpha=1; okPaso=false;
            while alpha>1/2048
                yy=y+alpha*delta;
                if all(yy>=-1e-10) && 1-p.K*yy(5)>0 && 1+yy(3)>0
                    Fy=modelo(yy,p); nf=nf+1;
                    if norm(yy-xn-h*Fy,inf)<norm(G,inf)
                        y=yy; okPaso=true; break;
                    end
                end
                alpha=alpha/2;
            end
            if ~okPaso, error('Newton: busqueda de paso fallo'); end
            itN=itN+1;
        end
        if ~ok
            F=modelo(y,p); nf=nf+1;
            if norm(y-xn-h*F,inf)>=o.tolG, error('Newton no convergio'); end
        end
        if any(~isfinite(y)) || any(y < -1e-10), error('Solucion no fisica'); end
        X(:,n+1)=y;
    end
end

function J=jac_f(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5);
    a=p.b2/(1-p.K*I); da=p.b2*p.K/(1-p.K*I)^2;
    b=p.beta*(1+M)/(1+Z); bm=p.beta/(1+Z); bz=-p.beta*(1+M)/(1+Z)^2; e=p.eps;
    J=[-p.b1*M-a*Z-p.mu1, -p.b1*U, -a*U, 0, -da*U*Z, 0;
       p.b1*M, p.b1*U-p.mu1, 0, 0, 0, 0;
       a*Z, 0, a*U-p.mu1, 0, da*U*Z, 0;
       0, -e*bm*S*I, -e*bz*S*I, e*(-b*I-p.mu2), -e*b*S, e*p.eta;
       0, e*bm*S*I, e*bz*S*I, e*b*I, e*(b*S-p.gamma-p.mu2), 0;
       0, 0, 0, 0, e*p.gamma, -e*(p.eta+p.mu2)];
end
