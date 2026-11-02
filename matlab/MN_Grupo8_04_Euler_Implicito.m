%% Grupo 8 - UMZ-SIRS - Euler hacia atras + Newton + Gauss-Seidel (VERSION CANONICA)
% CC2104 Metodos Numericos - UTEC. Ejecutar Run con matlab/ como directorio.
% Integrador: Euler hacia atras. Solver no lineal: Newton multivariable.
% Solver lineal interno: Gauss-Seidel (o solucion directa, solo como CONTROL).
% Esta es la implementacion unica; los archivos 05, 06 y 11 usan las mismas
% funciones locales (euler_atras, paso_implicito, gauss_seidel).
% Cada paso rechazado se clasifica por causa:
%   GS (rho>=1, residuo lineal o diagonal casi nula), positividad/conservacion,
%   Newton (no converge) y busqueda lineal (no reduce el residuo).
clear; clc; close all; format long g

%% 1. Parametros (iguales que en archivos 01, 02 y 03)
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=1e-2;
condiciones={[0.8;0.2;0;0.8;0.2;0], ...
             [0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original, Z(0)=0','Alternativa, Z(0)=0.1'};
T=100; h=0.1;
opt.tolResidual=1e-10;            % residuo no lineal ||G||_inf
opt.tolGS=1e-6;                   % cambio entre barridos de Gauss-Seidel
opt.tolLin=1e-11;                 % residuo lineal ||W y - b||_inf (se supervisa)
opt.maxNewton=50; opt.maxGS=50;
opt.tolPos=1e-10; opt.tolCons=1e-8;
opt.minH=1e-8; opt.maxRechazos=50;
opt.reducir=true;                 % true: rechazo y reduccion del paso
opt.solver='gs';                  % 'gs' o 'directo' (control)

%% 2. Integracion con h=0.1 en los dos estados iniciales
for caso=1:numel(condiciones)
    [t,X,info]=euler_atras(condiciones{caso},p,T,h,opt);
    ec=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
    mn=min(X,[],'all'); nn=nnz(X < -opt.tolPos);
    fprintf('\n===== Euler atras (Gauss-Seidel), h=%g: %s =====\n',h,nombres{caso});
    imprimir_info(info);
    fprintf('Error conservacion=%.3e, minimo=%.3e, negativos=%d\n',ec,mn,nn);
    fprintf('Estado final [U M Z S I R]:\n'); disp(X(:,end).');

    figure('Name',['Euler implicito - ' nombres{caso}],'Color','w');
    subplot(1,2,1);
    plot(t,X(1:3,:)','LineWidth',1.3); grid on;
    xlabel('Tiempo rápido t'); ylabel('Proporción'); title('Capa informativa (Euler atrás)');
    legend('U','M','Z','Location','best');
    subplot(1,2,2);
    plot(p.eps*t,X(4:6,:)','LineWidth',1.3); grid on;
    xlabel('Tiempo lento \tau=\epsilon t'); ylabel('Proporción'); title('Capa epidemiológica (Euler atrás)');
    legend('S','I','R','Location','best');
end

%% 3. Control: Gauss-Seidel frente a solucion directa con paso nominal grande
% Si con solucion directa desaparecen los rechazos, la limitacion es de
% Gauss-Seidel (rho(T_GS)>=1 o residuo lineal), no del integrador implicito.
hGrande=4; x0=condiciones{2};
solvers={'gs','directo'};
for k=1:2
    o=opt; o.solver=solvers{k};
    [t,X,info]=euler_atras(x0,p,T,hGrande,o); %#ok<ASGLU>
    fprintf('\n===== Euler atras, h nominal=%g, solver lineal: %s =====\n',hGrande,solvers{k});
    imprimir_info(info);
    if strcmp(info.estado,'OK')
        fprintf('h efectivo medio=%.4g (T/pasos aceptados)\n',T/info.aceptados);
    end
end
fprintf('\nEULER HACIA ATRAS COMPLETADO.\n');

%% FUNCIONES LOCALES
function imprimir_info(info)
    fprintf('Estado=%s\n',info.estado);
    fprintf('Pasos aceptados=%d, rechazados=%d (GS=%d, positividad=%d, Newton=%d, linea=%d)\n', ...
        info.aceptados,info.rechazados,info.rechGS,info.rechPos,info.rechNewton,info.rechLinea);
    fprintf('Newton=%d, barridos GS=%d, Jacobianas=%d, eval.f=%d\n', ...
        info.newton,info.gs,info.jac,info.evalf);
    fprintf('rho(T_GS) maximo observado=%.4f, h minimo aceptado=%.4g\n',info.rhoMax,info.hMin);
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
            % Condicion de convergencia de Gauss-Seidel: rho(T_GS)<1.
            Bgs=-(tril(W)\triu(W,1));
            rho=max(abs(eig(Bgs))); s.rho=max(s.rho,rho);
            if rho>=1, causa='gs_rho'; return, end
            [delta,exito,it]=gauss_seidel(W,-G,opt); s.gs=s.gs+it;
            if ~exito, causa='gs_residuo'; return, end
        else
            delta=W\(-G);
            if any(~isfinite(delta)), causa='singular'; return, end
        end
        % Amortiguacion: mantiene el dominio y exige reducir el residuo.
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
        % Se exige tambien el residuo lineal, no solo el tamano de la actualizacion.
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
