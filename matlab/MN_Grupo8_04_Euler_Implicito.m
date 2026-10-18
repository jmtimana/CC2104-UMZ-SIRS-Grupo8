%% Grupo 8 - UMZ-SIRS - Euler hacia atras + Newton + Gauss-Seidel
% CC2104 Metodos Numericos - UTEC
% Integrador: Euler hacia atras. Solver no lineal: Newton multivariable.
% Solver lineal interno: Gauss-Seidel. Ejecutar Run.
clear; clc; close all; format long g

%% 1. Parametros (iguales que en archivos 01, 02 y 03)
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=1e-2;
condiciones={[0.8;0.2;0;0.8;0.2;0], ...
             [0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original, Z(0)=0','Alternativa, Z(0)=0.1'};
T=100; h=0.1;
% Tolerancias iniciales del proyecto + control adicional de residuos.
opt.tolNewton=1e-8; opt.tolResidual=1e-10;
opt.tolGS=1e-6; opt.tolResiduoLineal=1e-11;
opt.maxNewton=50; opt.maxGS=50;
opt.tolPos=1e-10; opt.tolCons=1e-8;
opt.minH=1e-8; opt.maxRechazos=30;

%% 2. Integracion y metricas
for caso=1:numel(condiciones)
    [t,X,info]=integrar_atras(condiciones{caso},p,T,h,opt);
    e1=max(abs(sum(X(1:3,:),1)-1));
    e2=max(abs(sum(X(4:6,:),1)-1));
    ec=max(e1,e2); mn=min(X,[],'all');
    nn=nnz(X < -opt.tolPos);
    fprintf('\n===== Euler atras: %s =====\n',nombres{caso});
    fprintf('Pasos aceptados=%d, rechazados=%d\n', ...
        info.aceptados,info.rechazados);
    fprintf('Iteraciones Newton=%d, iteraciones GS=%d\n', ...
        info.iterNewton,info.iterGS);
    fprintf('Max residuo no lineal=%.3e\n',info.maxResidual);
    fprintf('Error conservacion=%.3e, minimo=%.3e, negativos=%d\n',ec,mn,nn);
    fprintf('Conservacion: %s | Positividad: %s\n', ...
        estado(ec<opt.tolCons),estado(nn==0));
    fprintf('Estado final [U M Z S I R]:\n'); disp(X(:,end).');

    figure('Name',['Euler implicito - ' nombres{caso}],'Color','w');
    subplot(1,2,1);
    plot(t,X(1:3,:)','LineWidth',1.3); grid on;
    xlabel('t rapido'); ylabel('Proporcion'); title('Informacion: Euler atras');
    legend('U','M','Z','Location','best');
    subplot(1,2,2);
    plot(p.eps*t,X(4:6,:)','LineWidth',1.3); grid on;
    xlabel('tau = epsilon*t'); ylabel('Proporcion'); title('Epidemia: Euler atras');
    legend('S','I','R','Location','best');
end
fprintf('\nEULER HACIA ATRAS COMPLETADO.\n');

%% FUNCIONES LOCALES
function [t,X,info]=integrar_atras(x0,p,T,hbase,opt)
    t=0; X=x0;
    info.aceptados=0; info.rechazados=0;
    info.iterNewton=0; info.iterGS=0; info.maxResidual=0;
    h=hbase;
    while t(end)<T
        h=min(h,T-t(end));
        xAnterior=X(:,end);
        [xNuevo,ok,nNewton,nGS,r]=paso_implicito(xAnterior,p,h,opt);
        info.iterNewton=info.iterNewton+nNewton;
        info.iterGS=info.iterGS+nGS;
        if ok && all(xNuevo>=-opt.tolPos) && ...
                abs(sum(xNuevo(1:3))-1)<opt.tolCons && ...
                abs(sum(xNuevo(4:6))-1)<opt.tolCons
            t(end+1)=t(end)+h; %#ok<AGROW>
            X(:,end+1)=xNuevo; %#ok<AGROW>
            info.aceptados=info.aceptados+1;
            info.maxResidual=max(info.maxResidual,r);
            h=min(hbase,1.5*h); % recuperar paso tras rechazos
        else
            info.rechazados=info.rechazados+1;
            h=h/2;
            if h<opt.minH || info.rechazados>opt.maxRechazos
                error('Euler atras: no se pudo resolver un paso: revisar GS/Newton.');
            end
        end
    end
end

function [x,ok,nNewton,nGS,r]=paso_implicito(xn,p,h,opt)
    x=xn; nGS=0; ok=false; r=inf; nNewton=0;
    for k=1:opt.maxNewton
        nNewton=k;
        G=x-xn-h*f(x,p);
        r=norm(G,inf);
        if r<opt.tolResidual
            ok=true; return
        end
        W=eye(6)-h*jacobiana_f(x,p);
        if any(~isfinite(W),'all') || any(abs(diag(W))<1e-13)
            return
        end
        % Comprobar la convergencia de Gauss-Seidel.
        Bgs=-(tril(W)\triu(W,1));
        if max(abs(eig(Bgs)))>=1
            return
        end
        [delta,exito,itGS]=gauss_seidel(W,-G,opt);
        nGS=nGS+itGS;
        if ~exito, return; end
        % Amortiguacion para mantener el dominio y reducir el residual.
        alpha=1; aceptado=false;
        while alpha>=1/1024
            candidato=x+alpha*delta;
            if all(isfinite(candidato)) && all(candidato>=-opt.tolPos) ...
                    && 1-p.K*candidato(5)>0 && 1+candidato(3)>0
                nuevoG=candidato-xn-h*f(candidato,p);
                if norm(nuevoG,inf)<r
                    x=candidato; aceptado=true; break;
                end
            end
            alpha=alpha/2;
        end
        if ~aceptado, return; end
        if norm(alpha*delta,inf)<opt.tolNewton && ...
                norm(x-xn-h*f(x,p),inf)<opt.tolResidual
            ok=true; r=norm(x-xn-h*f(x,p),inf); return
        end
    end
    r=norm(x-xn-h*f(x,p),inf);
end

function [y,ok,k]=gauss_seidel(A,b,opt)
    n=numel(b); y=zeros(n,1); ok=false;
    for k=1:opt.maxGS
        ant=y;
        for i=1:n
            y(i)=(b(i)-A(i,1:i-1)*y(1:i-1) ...
                    -A(i,i+1:end)*ant(i+1:end))/A(i,i);
        end
        res=norm(A*y-b,inf);
        % La tolerancia antigua (1e-6) controla cambios; el residuo
        % estricto evita detener GS antes de que Newton lo necesite.
        if norm(y-ant,inf)<opt.tolGS && ...
                res<opt.tolResiduoLineal
            ok=true; return
        end
    end
end

function J=jacobiana_f(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5);
    a=p.b2/(1-p.K*I);
    da=p.b2*p.K/(1-p.K*I)^2;
    b=p.beta*(1+M)/(1+Z);
    bm=p.beta/(1+Z);
    bz=-p.beta*(1+M)/(1+Z)^2;
    e=p.eps;
    J=[-p.b1*M-a*Z-p.mu1, -p.b1*U, -a*U, 0, -da*U*Z, 0;
        p.b1*M, p.b1*U-p.mu1, 0, 0, 0, 0;
        a*Z, 0, a*U-p.mu1, 0, da*U*Z, 0;
        0, -e*bm*S*I, -e*bz*S*I, e*(-b*I-p.mu2), -e*b*S, e*p.eta;
        0, e*bm*S*I, e*bz*S*I, e*b*I, e*(b*S-p.gamma-p.mu2), 0;
        0, 0, 0, 0, e*p.gamma, -e*(p.eta+p.mu2)];
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
