%% GRUPO 8 - Comparacion cuantitativa de integradores UMZ-SIRS (ENTREGA FINAL)
% MATLAB R2019b+ | Ejecutar Run con matlab/ como directorio de trabajo.
% Solucion de referencia: ode45 con tolerancias estrictas (NO exacta; se valida
% en el archivo 09). Pasos fijos: Euler explicito, RK4, Euler atras (canonico,
% sin reduccion de paso). RKF45: se barre la TOLERANCIA (no el paso inicial).
% Salidas: ../resultados/MN_Grupo8_resultados_caso_1.csv y _caso_2.csv
%          ../figuras/fig01_convergencia_alternativa.png
%          ../figuras/fig10_trabajo_precision.png
% Nota: guardar en UTF-8 (los textos de las figuras llevan tildes).
clear; clc; close all; format long g

carpetaFiguras='../figuras'; carpetaRes='../resultados';

% LETRA LEGIBLE AL ANCHO DE COLUMNA DEL PAPER

set(groot,'defaultAxesFontSize',14,'defaultLegendFontSize',12,'defaultTextFontSize',14);

if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Parametros y controles
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=1e-2;
T=100;
hs=[0.8 0.4 0.2 0.1 0.05 0.025];      % todos dividen exactamente T
rtols=10.^(-4:-1:-10);                % barrido de tolerancias de RKF45
condiciones={[0.8;0.2;0;0.8;0.2;0],[0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original Z0=0','Alternativa Z0=0.1'};
metodosFijos={'Euler explicito','RK4','Euler atras'};
optRef=odeset('RelTol',1e-12,'AbsTol',1e-13,'MaxStep',0.01);

opt.tolResidual=1e-10; opt.tolGS=1e-6; opt.tolLin=1e-11;
opt.maxNewton=50; opt.maxGS=50; opt.tolPos=1e-10; opt.tolCons=1e-8;
opt.minH=1e-8; opt.maxRechazos=50;
opt.reducir=false;                    % paso fijo: si falla, se registra FALLO
opt.solver='gs';
opt.hmin=1e-10; opt.maxIntentos=200000;
h0RKF=0.1; hmaxRKF=2;
nRep=5;                               % repeticiones para medir tiempo (se usa la mediana)

%% 2. Comparacion en cada condicion inicial
for caso=1:2
    x0=condiciones{caso};
    fprintf('\n================ %s ================\n',nombres{caso});
    reloj=tic; solRef=ode45(@(t,x)modelo(x,p),[0 T],x0,optRef);
    fprintf('Referencia ode45: %.3f s, %d nodos internos\n',toc(reloj),numel(solRef.x));
    filas={};
    % --- Integradores de paso fijo
    for im=1:3
        errPrev=NaN;
        for ih=1:numel(hs)
            h=hs(ih);
            segs=zeros(1,nRep);
            for rep=1:nRep+1                  % la 1.a corrida es de calentamiento y se descarta
                reloj=tic;
                switch im
                    case 1, [t,X,info]=explicito(x0,p,T,h,1,opt);
                    case 2, [t,X,info]=explicito(x0,p,T,h,4,opt);
                    case 3, [t,X,info]=euler_atras(x0,p,T,h,opt);
                end
                if rep>1, segs(rep-1)=toc(reloj); end
            end
            seg=median(segs); info=completar(info);
            [err,cons,mn]=metricas(t,X,solRef,T,info);
            orden=NaN;
            if ih>1 && isfinite(err) && isfinite(errPrev) && err>0
                orden=log(errPrev/err)/log(hs(ih-1)/h);   % orden observado
            end
            errPrev=err;
            filas(end+1,:)={caso,metodosFijos{im},h,NaN,err,orden,cons,mn, ...
                info.aceptados,info.rechazados,info.evalf,info.jac,info.newton, ...
                info.gs,seg,info.estado}; %#ok<AGROW>
        end
    end
    % --- RKF45: barrido de tolerancia
    for k=1:numel(rtols)
        rt=rtols(k);
        segs=zeros(1,nRep);
        for rep=1:nRep+1                      % la 1.a corrida es de calentamiento y se descarta
            reloj=tic;
            [t,X,info]=rkf45(x0,p,T,h0RKF,hmaxRKF,rt,1e-2*rt,opt);
            if rep>1, segs(rep-1)=toc(reloj); end
        end
        seg=median(segs); info=completar(info);
        [err,cons,mn]=metricas(t,X,solRef,T,info);
        filas(end+1,:)={caso,'RKF45',NaN,rt,err,NaN,cons,mn,info.aceptados, ...
            info.rechazados,info.evalf,0,0,0,seg,info.estado}; %#ok<AGROW>
    end
    tabla=cell2table(filas,'VariableNames', ...
        {'Caso','Metodo','h','RelTol','ErrorMaxRef','OrdenObs','ErrorConservacion', ...
         'Minimo','Pasos','Rechazados','EvalF','Jacobianas','IterNewton','IterGS', ...
         'Segundos','Estado'});
    disp(tabla);
    nombreCSV=sprintf('MN_Grupo8_resultados_caso_%d.csv',caso);
    writetable(tabla,fullfile(carpetaRes,nombreCSV));
    fprintf('Tabla guardada: %s\n',nombreCSV);

    if caso==2
        %% Figura 01: convergencia (paso fijo) con pendientes de referencia
        fig1=figure('Color','w','Name','Figura 01 - Convergencia','Position',[60 60 960 510]);
        hold on;
        estilos={'o-','s--','^-.'}; colores=lines(3); hl=gobjects(1,3);
        for im=1:3
            m=strcmp(tabla.Metodo,metodosFijos{im}) & strncmp(tabla.Estado,'OK',2);
            hl(im)=plot(tabla.h(m),max(tabla.ErrorMaxRef(m),eps),estilos{im}, ...
                'Color',colores(im,:),'LineWidth',1.5,'MarkerSize',7);
        end
        ie=strcmp(tabla.Metodo,'Euler explicito') & tabla.h==0.2;
        ir=strcmp(tabla.Metodo,'RK4') & tabla.h==0.2;
        hh=[min(hs) max(hs)];
        r1=plot(hh,tabla.ErrorMaxRef(ie)*(hh/0.2),'k:','LineWidth',1.2);
        r4=plot(hh,tabla.ErrorMaxRef(ir)*(hh/0.2).^4,'k--','LineWidth',1.2);
        set(gca,'XScale','log','YScale','log'); grid on;
        xticks(sort(hs)); xticklabels(string(sort(hs))); set(gca,'XMinorTick','off');
        xlabel('Paso fijo h'); ylabel('Error máximo frente a la referencia');
        yline(1e-11,'-.','incertidumbre de la referencia','Color',[0.4 0.4 0.4], ...
            'LineWidth',1.1,'LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
        title(sprintf('Convergencia de los integradores de paso fijo (\\epsilon=%g, Z(0)=0.1)',p.eps));
        legend([hl r1 r4],{'Euler explícito','RK4','Euler atrás','Pendiente O(h)', ...
            'Pendiente O(h^4)'},'Location','southeast');
        exportgraphics(fig1,fullfile(carpetaFiguras,'fig01_convergencia_alternativa.png'),'Resolution',300);
        fprintf('Figura exportada: fig01_convergencia_alternativa.png\n');

        %% Figura 10: trabajo-precision (error vs evaluaciones de f y vs tiempo)
        fig10=figure('Color','w','Name','Figura 10 - Trabajo-precision','Position',[60 60 1150 450]);
        metodosTP={'Euler explicito','RK4','Euler atras','RKF45'};
        etiq={'Euler explícito','RK4','Euler atrás (Gauss–Seidel)','RKF45 (tol. variable)'};
        col=lines(4); mk={'o-','s-','^-','d-'};
        for panel=1:2
            subplot(1,2,panel); hold on;
            for im=1:4
                m=strcmp(tabla.Metodo,metodosTP{im}) & strncmp(tabla.Estado,'OK',2);
                if panel==1, xx=tabla.EvalF(m); else, xx=tabla.Segundos(m); end
                plot(xx,max(tabla.ErrorMaxRef(m),eps),mk{im},'Color',col(im,:), ...
                    'LineWidth',1.4,'MarkerSize',6);
            end
            set(gca,'XScale','log','YScale','log'); grid on;
            ylabel('Error máximo frente a la referencia');
            if panel==1
                xlabel('Evaluaciones de f'); title('(a) Error frente a evaluaciones de f');
                legend(etiq,'Location','northeast');
            else
                xlabel('Tiempo de cómputo (s)'); title('(b) Error frente al tiempo');
            end
        end
        exportgraphics(fig10,fullfile(carpetaFiguras,'fig10_trabajo_precision.png'),'Resolution',300);
        fprintf('Figura exportada: fig10_trabajo_precision.png\n');
    end
end
fprintf('\nCOMPARACION COMPLETADA. Revisar posibles FALLOS y discrepancias.\n');
fprintf('El tiempo de ejecucion depende de MATLAB, equipo y carga de trabajo.\n');
fprintf('EvalF cuenta solo evaluaciones de f; Euler atras ademas usa Jacobianas y barridos GS.\n');

%% FUNCIONES LOCALES
function [err,cons,mn]=metricas(t,X,solRef,T,info)
    if strcmp(info.estado,'OK')
        err=max(abs(X-deval(solRef,min(max(t,0),T))),[],'all');
        cons=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
        mn=min(X,[],'all');
    else
        err=NaN; cons=NaN; mn=NaN;
    end
end

function info=completar(info)
    campos={'jac','newton','gs','rechGS','rechPos','rechNewton','rechLinea','rhoMax','evalf'};
    for k=1:numel(campos)
        if ~isfield(info,campos{k}), info.(campos{k})=0; end
    end
end

function [t,X,info]=explicito(x0,p,T,hBase,orden,opt)
    cap=1024; tb=zeros(1,cap); Xb=zeros(6,cap); n=1; Xb(:,1)=x0;
    h=hBase; seguidos=0;
    info=struct('estado','OK','aceptados',0,'rechazados',0,'evalf',0,'rechPos',0,'hMin',inf);
    while T-tb(n)>1e-10*max(1,T)
        h=min(h,T-tb(n)); x=Xb(:,n);
        if orden==1
            xn=x+h*modelo(x,p); info.evalf=info.evalf+1;
        else
            k1=modelo(x,p); k2=modelo(x+h*k1/2,p);
            k3=modelo(x+h*k2/2,p); k4=modelo(x+h*k3,p);
            xn=x+h*(k1+2*k2+2*k3+k4)/6; info.evalf=info.evalf+4;
        end
        if ~all(isfinite(xn)) || any(abs(xn)>1e8), causa='no finito';
        elseif any(xn<-opt.tolPos), causa='negativo';
        elseif abs(sum(xn(1:3))-1)>=opt.tolCons || abs(sum(xn(4:6))-1)>=opt.tolCons, causa='conservacion';
        elseif 1-p.K*xn(5)<=0 || 1+xn(3)<=0, causa='dominio';
        else, causa='ok';
        end
        if strcmp(causa,'ok')
            n=n+1; if n>cap, cap=2*cap; tb(cap)=0; Xb(6,cap)=0; end
            tb(n)=tb(n-1)+h; Xb(:,n)=xn;
            info.aceptados=info.aceptados+1; info.hMin=min(info.hMin,h);
            h=min(hBase,1.5*h); seguidos=0;
        else
            info.rechazados=info.rechazados+1; info.rechPos=info.rechPos+1; seguidos=seguidos+1;
            if ~opt.reducir
                info.estado=['INADMISIBLE (' causa ')']; break
            end
            h=h/2;
            if h<opt.minH || seguidos>opt.maxRechazos
                info.estado=['INADMISIBLE (' causa ')']; break
            end
        end
    end
    t=tb(1:n); X=Xb(:,1:n);
end

function [t,X,info]=rkf45(x0,p,T,h0,hmax,rtol,atol,opt)
    cap=1024; tb=zeros(1,cap); Xb=zeros(6,cap); n=1; Xb(:,1)=x0;
    h=h0; intentos=0;
    info=struct('estado','OK','aceptados',0,'rechazados',0,'evalf',0,'rechPos',0,'hMin',inf,'hMax',0);
    while T-tb(n)>1e-12*max(1,T)
        intentos=intentos+1;
        if intentos>opt.maxIntentos || h<opt.hmin
            info.estado='FALLO (paso minimo o demasiados intentos)'; break
        end
        h=min([h,hmax,T-tb(n)]); x=Xb(:,n);
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
        er=max(abs(x5-x4)./esc); if ~isfinite(er), er=inf; end
        valido=all(isfinite(x5)) && all(x5>=-opt.tolPos) && 1-p.K*x5(5)>0 && 1+x5(3)>0;
        if er<=1 && valido
            n=n+1; if n>cap, cap=2*cap; tb(cap)=0; Xb(6,cap)=0; end
            tb(n)=tb(n-1)+h; Xb(:,n)=x5;
            info.aceptados=info.aceptados+1;
            info.hMin=min(info.hMin,h); info.hMax=max(info.hMax,h);
        else
            info.rechazados=info.rechazados+1;
            if ~valido, info.rechPos=info.rechPos+1; end
        end
        if er==0, fac=2; else, fac=min(2,max(0.2,0.9*(1/er)^(1/5))); end
        if ~valido, fac=min(fac,0.5); end
        h=h*fac;
    end
    t=tb(1:n); X=Xb(:,1:n);
end

function [t,X,info]=euler_atras(x0,p,T,hBase,opt)
    cap=1024; tb=zeros(1,cap); Xb=zeros(6,cap); n=1; Xb(:,1)=x0;
    h=hBase; seguidos=0;
    info=struct('estado','OK','aceptados',0,'rechazados',0,'newton',0,'gs',0, ...
        'jac',0,'evalf',0,'rechGS',0,'rechPos',0,'rechNewton',0,'rechLinea',0, ...
        'rhoMax',0,'hMin',inf);
    while T-tb(n)>1e-10*max(1,T)
        h=min(h,T-tb(n));
        [xn,ok,causa,s]=paso_implicito(Xb(:,n),p,h,opt);
        info.newton=info.newton+s.newton; info.gs=info.gs+s.gs;
        info.jac=info.jac+s.jac; info.evalf=info.evalf+s.evalf;
        info.rhoMax=max(info.rhoMax,s.rho);
        if ok && ~(all(xn>=-opt.tolPos) && abs(sum(xn(1:3))-1)<opt.tolCons ...
                && abs(sum(xn(4:6))-1)<opt.tolCons)
            ok=false; causa='positividad';
        end
        if ok
            n=n+1; if n>cap, cap=2*cap; tb(cap)=0; Xb(6,cap)=0; end
            tb(n)=tb(n-1)+h; Xb(:,n)=xn;
            info.aceptados=info.aceptados+1; info.hMin=min(info.hMin,h);
            h=min(hBase,1.5*h); seguidos=0;
        else
            info.rechazados=info.rechazados+1; seguidos=seguidos+1;
            switch causa
                case {'gs_rho','gs_residuo'}, info.rechGS=info.rechGS+1;
                case 'positividad',           info.rechPos=info.rechPos+1;
                case 'linea',                 info.rechLinea=info.rechLinea+1;
                otherwise,                    info.rechNewton=info.rechNewton+1;
            end
            if ~opt.reducir
                info.estado=['FALLO (' causa ')']; break
            end
            h=h/2;
            if h<opt.minH || seguidos>opt.maxRechazos
                info.estado=['FALLO (' causa ')']; break
            end
        end
    end
    t=tb(1:n); X=Xb(:,1:n);
end

function [x,ok,causa,s]=paso_implicito(xn,p,h,opt)
    x=xn; ok=false; causa='newton_maxiter';
    s=struct('newton',0,'gs',0,'jac',0,'evalf',0,'rho',0);
    for k=1:opt.maxNewton
        G=x-xn-h*modelo(x,p); s.evalf=s.evalf+1;
        r=norm(G,inf);
        if r<opt.tolResidual, ok=true; causa='ok'; return, end
        W=eye(6)-h*jacobiana_f(x,p); s.jac=s.jac+1;
        if any(~isfinite(W),'all') || any(abs(diag(W))<1e-13)
            causa='singular'; return
        end
        if strcmp(opt.solver,'gs')
            Bgs=-(tril(W)\triu(W,1));
            rho=max(abs(eig(Bgs))); s.rho=max(s.rho,rho);
            if rho>=1, causa='gs_rho'; return, end
            [delta,exito,it]=gauss_seidel(W,-G,opt); s.gs=s.gs+it;
            if ~exito, causa='gs_residuo'; return, end
        else
            delta=W\(-G);
            if any(~isfinite(delta)), causa='singular'; return, end
        end
        alpha=1; aceptado=false;
        while alpha>=1/1024
            cand=x+alpha*delta;
            if all(isfinite(cand)) && all(cand>=-opt.tolPos) && ...
                    1-p.K*cand(5)>0 && 1+cand(3)>0
                v=cand-xn-h*modelo(cand,p); s.evalf=s.evalf+1;
                if norm(v,inf)<r, x=cand; aceptado=true; break, end
            end
            alpha=alpha/2;
        end
        if ~aceptado, causa='linea'; return, end
        s.newton=s.newton+1;
    end
    G=x-xn-h*modelo(x,p); s.evalf=s.evalf+1;
    if norm(G,inf)<opt.tolResidual, ok=true; causa='ok'; end
end

function [y,ok,k]=gauss_seidel(A,b,opt)
    n=numel(b); y=zeros(n,1); ok=false;
    for k=1:opt.maxGS
        ant=y;
        for i=1:n
            y(i)=(b(i)-A(i,1:i-1)*y(1:i-1)-A(i,i+1:n)*ant(i+1:n))/A(i,i);
        end
        if norm(y-ant,inf)<opt.tolGS && norm(A*y-b,inf)<opt.tolLin
            ok=true; return
        end
    end
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

function J=jacobiana_f(x,p)
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
