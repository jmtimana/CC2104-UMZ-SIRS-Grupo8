%% Grupo 8 - Integracion temporal del modelo UMZ-SIRS: Euler y RK4
% Metodos Numericos CC2104 - UTEC, 2026-2.
% Ejecutar este script completo (Run); tambien se puede abrir en Live Editor.
% Se conserva el mismo modelo y parametros del archivo 01.
clear; clc; close all;
format long g

%% 1. Parametros exploratorios (NO atribuir al articulo original)
p.b1=1.1; p.b2=0.92; p.mu1=1.0; p.K=0.9;
p.beta=6.0; p.gamma=1.2; p.eta=0.7; p.mu2=1.0; p.eps=1e-2;
variables={'U','M','Z','S','I','R'};

%% 2. Configuracion de los experimentos
% t es la escala rapida; tau=eps*t es la escala epidemiologica.
T=100;                    % horizonte en la escala rapida (tau_final=1)
h=0.1;                    % paso temporal de ambos metodos
assert(h>0 && T>0 && abs(T/h-round(T/h))<1e-10, ...
    'T/h debe ser entero para comparar las mismas mallas.');
condiciones={[0.8;0.2;0;0.8;0.2;0], ...
             [0.7;0.2;0.1;0.8;0.2;0]};
nombres={'Original (Z(0)=0)','Alternativa (Z(0)=0.1)'};
tol_cons=1e-8;
tol_pos=1e-10;

%% 3. Integracion y verificacion en los dos escenarios iniciales
for caso=1:numel(condiciones)
    x0=condiciones{caso};
    assert(all(x0>=0) && abs(sum(x0(1:3))-1)<1e-12 && ...
        abs(sum(x0(4:6))-1)<1e-12, 'Condicion inicial invalida.');

    [t,Xe,evalEuler]=integrar_umz(@paso_euler,x0,p,T,h);
    [~,Xr,evalRK4]=integrar_umz(@paso_rk4,x0,p,T,h);
    me=metricas(Xe,tol_cons,tol_pos);
    mr=metricas(Xr,tol_cons,tol_pos);

    % Comparacion de refinamiento h vs h/2. No es error exacto.
    [~,Xe2]=integrar_umz(@paso_euler,x0,p,T,h/2);
    [~,Xr2]=integrar_umz(@paso_rk4,x0,p,T,h/2);
    difEuler=max(abs(Xe-Xe2(:,1:2:end)),[],'all');
    difRK4=max(abs(Xr-Xr2(:,1:2:end)),[],'all');

    fprintf('\n======= %s =======\n',nombres{caso});
    fprintf('T=%.3g, h=%.3g, pasos=%d, tau_final=%.3g\n',T,h,numel(t)-1,p.eps*T);
    fprintf('Metodo          Eval.f     Error cons.      Minimo       Negativos    Dif.h-h/2\n');
    fprintf('Euler       %10d    %.3e      % .3e    %7d      %.3e\n', ...
        evalEuler,me.errCons,me.minimo,me.negativos,difEuler);
    fprintf('RK4         %10d    %.3e      % .3e    %7d      %.3e\n', ...
        evalRK4,mr.errCons,mr.minimo,mr.negativos,difRK4);
    fprintf('Conservacion Euler: %s | RK4: %s\n',estado(me.consOK),estado(mr.consOK));
    fprintf('Positividad Euler: %s | RK4: %s\n',estado(me.posOK),estado(mr.posOK));

    %% 4. Graficas de informacion y epidemia
    figure('Name',['UMZ-SIRS - ' nombres{caso}],'Color','w');
    subplot(2,2,1);
    plot(t,Xe(1:3,:)','LineWidth',1.25); grid on;
    title('Informacion: Euler explicito'); xlabel('t (rapido)'); ylabel('Proporcion');
    legend(variables(1:3),'Location','best');

    subplot(2,2,2);
    plot(t,Xr(1:3,:)','LineWidth',1.25); grid on;
    title('Informacion: RK4'); xlabel('t (rapido)'); ylabel('Proporcion');
    legend(variables(1:3),'Location','best');

    subplot(2,2,3);
    plot(p.eps*t,Xe(4:6,:)','LineWidth',1.25); grid on;
    title('Epidemia: Euler explicito'); xlabel('tau = epsilon*t (lento)'); ylabel('Proporcion');
    legend(variables(4:6),'Location','best');

    subplot(2,2,4);
    plot(p.eps*t,Xr(4:6,:)','LineWidth',1.25); grid on;
    title('Epidemia: RK4'); xlabel('tau = epsilon*t (lento)'); ylabel('Proporcion');
    legend(variables(4:6),'Location','best');

    %% 5. Comparacion de la poblacion infectada
    figure('Name',['Comparacion de I - ' nombres{caso}],'Color','w');
    plot(p.eps*t,Xe(5,:),'--','LineWidth',1.4); hold on;
    plot(p.eps*t,Xr(5,:),'-','LineWidth',1.4); grid on;
    xlabel('tau = epsilon*t'); ylabel('Fraccion infectada I');
    title(['Comparacion de Euler y RK4: ' nombres{caso}]);
    legend('Euler explicito','RK4','Location','best');
end

fprintf('\nINTEGRACION EULER/RK4 COMPLETADA.\n');
fprintf('Las diferencias h versus h/2 son estimaciones de refinamiento, no errores exactos.\n');

%% FUNCIONES LOCALES
function [t,X,evals]=integrar_umz(paso,x0,p,T,h)
    n=round(T/h);
    t=linspace(0,T,n+1);
    X=zeros(6,n+1); X(:,1)=x0;
    evals=0;
    for k=1:n
        [X(:,k+1),nf]=paso(X(:,k),h,p);
        evals=evals+nf;
        if any(~isfinite(X(:,k+1)))
            error('Valores no finitos en paso %d (t=%g).',k,t(k+1));
        end
    end
end

function [xnext,nf]=paso_euler(x,h,p)
    xnext=x+h*modelo_umz_sirs(x,p);
    nf=1;
end

function [xnext,nf]=paso_rk4(x,h,p)
    k1=modelo_umz_sirs(x,p);
    k2=modelo_umz_sirs(x+(h/2)*k1,p);
    k3=modelo_umz_sirs(x+(h/2)*k2,p);
    k4=modelo_umz_sirs(x+h*k3,p);
    xnext=x+(h/6)*(k1+2*k2+2*k3+k4);
    nf=4;
end

function m=metricas(X,tolCons,tolPos)
    e1=abs(sum(X(1:3,:),1)-1);
    e2=abs(sum(X(4:6,:),1)-1);
    m.errCons=max([e1,e2]);
    m.minimo=min(X,[],'all');
    m.negativos=nnz(X < -tolPos);
    m.consOK=m.errCons<tolCons;
    m.posOK=m.negativos==0;
end

function s=estado(ok)
    if ok, s='OK'; else, s='NO CUMPLE'; end
end

function dx=modelo_umz_sirs(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5); R=x(6);
    if 1-p.K*I<=0 || 1+Z<=0
        error('El modelo encontro un denominador no positivo. Reducir h o revisar parametros.');
    end
    b2bar=p.b2/(1-p.K*I);
    betabar=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-b2bar*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        b2bar*U*Z-p.mu1*Z;
        p.eps*(p.mu2-betabar*S*I+p.eta*R-p.mu2*S);
        p.eps*(betabar*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end
