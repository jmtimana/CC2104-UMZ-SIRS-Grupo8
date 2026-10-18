%% Grupo 8 - UMZ-SIRS | 07. Escenarios epidemiologicos
% Metodo numerico principal: RK4 implementado por el grupo.
% Validacion independiente: ode15s (NO se presenta como metodo propio).
% Escenarios I y II: parametros consignados en el borrador de proyecto.
% Condiciones iniciales y tau_final: decisiones exploratorias del grupo;
% NO se atribuyen al articulo de Bulai et al. (2024).
clear; clc; close all; format long g

%% Exportacion automatica para el paper
carpetaFiguras = 'figuras';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end

%% 1. Escenarios y condiciones iniciales
p.mu1=5.2143; p.mu2=5.2143; p.gamma=0.1667;
p.eta=0.1; p.K=0.9; p.eps=1/50;

tauFinal=30; T=tauFinal/p.eps;
h=0.1;                     % paso fijo en tiempo rapido
I0=0.02; S0=1-I0; R0=0;

% Escenario I: beta / (gamma+mu2)<1, rama desinformada >1
casos(1).nombre='I - Desinformacion';
casos(1).b1=6.7041; casos(1).b2=5.9041; casos(1).beta=4.8429;
casos(1).informacion=[0.8;0.2;0];
% Escenario II: sin informacion >1, rama esceptica <1
casos(2).nombre='II - Escepticismo';
casos(2).b1=5.9041; casos(2).b2=6.7041; casos(2).beta=5.919;
casos(2).informacion=[0.8;0;0.2];

% El caso de control U=1, M=Z=0 es un estado invariante exacto.
infoControl=[1;0;0];

%% 2. Integracion y comparacion con control
filas={};
for c=1:numel(casos)
    q=p; q.b1=casos(c).b1; q.b2=casos(c).b2; q.beta=casos(c).beta;
    R00=q.beta/(q.gamma+q.mu2);
    R01=R00*(2*q.b1-q.mu1)/q.b1;
    R02=R00*q.b2/(2*q.b2-q.mu1);
    fprintf('\n===== %s =====\n',casos(c).nombre);
    fprintf('Umbrales reportados en el borrador: R00=%.6f, R01=%.6f, R02=%.6f\n',R00,R01,R02);
    if c==1
        fprintf('Hipotesis: R00<1<R01: %s\n',si_no(R00<1 && R01>1));
    else
        fprintf('Hipotesis: R02<1<R00: %s\n',si_no(R02<1 && R00>1));
    end

    x0control=[infoControl;S0;I0;R0];
    x0acoplado=[casos(c).informacion;S0;I0;R0];
    [tC,XC,mC]=rk4_integrar(x0control,q,T,h);
    [tA,XA,mA]=rk4_integrar(x0acoplado,q,T,h);
    if ~mC.valido || ~mA.valido
        warning('Caso %d: RK4 dejo dominio admisible. Se omite conclusion fisica.',c);
        continue
    end
    % Se coteja el error global del RK4 contra ode15s en sus mismos nodos.
    opciones=odeset('RelTol',1e-11,'AbsTol',1e-13,'MaxStep',1);
    solC=ode15s(@(t,x)f_umzsirs(x,q),[0 T],x0control,opciones);
    solA=ode15s(@(t,x)f_umzsirs(x,q),[0 T],x0acoplado,opciones);
    errC=max(abs(XC-deval(solC,tC)),[],'all');
    errA=max(abs(XA-deval(solA,tA)),[],'all');
    % Promedio en el ultimo 10 %% del tiempo lento: no confundir con limite t->infinito.
    ultC=tC>=0.9*T; ultA=tA>=0.9*T;
    IfinC=XC(5,end); IfinA=XA(5,end);
    IpromC=mean(XC(5,ultC)); IpromA=mean(XA(5,ultA));
    picoC=max(XC(5,:)); picoA=max(XA(5,:));
    fprintf('RK4, h=%g, tau final=%g: %d pasos por trayectoria\n',h,tauFinal,numel(tC)-1);
    fprintf('Control: I_final=%.8f, I_prom_ultimo10%%=%.8f, I_pico=%.8f\n',IfinC,IpromC,picoC);
    fprintf('Informacion: I_final=%.8f, I_prom_ultimo10%%=%.8f, I_pico=%.8f\n',IfinA,IpromA,picoA);
    fprintf('Diferencia I_final (informacion-control) = %.8f\n',IfinA-IfinC);
    fprintf('Error max RK4 vs ode15s: control %.3e | informacion %.3e\n',errC,errA);
    fprintf('Conservacion max: control %.3e | informacion %.3e\n',mC.cons,mA.cons);
    fprintf('Minima componente: control %.3e | informacion %.3e\n',mC.minimo,mA.minimo);

    filas(end+1,:)={c,'Control',R00,R01,R02,IfinC,IpromC,picoC,errC,mC.cons,mC.minimo}; %#ok<AGROW>
    filas(end+1,:)={c,'Con informacion',R00,R01,R02,IfinA,IpromA,picoA,errA,mA.cons,mA.minimo}; %#ok<AGROW>

    figure('Color','w','Name',casos(c).nombre);
    subplot(2,2,1)
    plot(q.eps*tC,XC(5,:),'--','LineWidth',1.4); hold on
    plot(q.eps*tA,XA(5,:),'LineWidth',1.4); grid on
    xlabel('Tiempo lento tau'); ylabel('Fraccion infectada I')
    title('Contraste epidemiologico'); legend('Sin informacion','Con informacion','Location','best')
    subplot(2,2,2)
    plot(q.eps*tA,XA(1:3,:)','LineWidth',1.2); grid on
    xlabel('Tiempo lento tau'); ylabel('Proporcion')
    title('Capa informativa del escenario'); legend('U','M','Z','Location','best')
    subplot(2,2,3)
    plot(q.eps*tA,XA(4:6,:)','LineWidth',1.2); grid on
    xlabel('Tiempo lento tau'); ylabel('Proporcion')
    title('SIRS con informacion'); legend('S','I','R','Location','best')
    subplot(2,2,4)
    plot(q.eps*tA,XA(5,:)-XC(5,:),'LineWidth',1.3); yline(0,'--'); grid on
    xlabel('Tiempo lento tau'); ylabel('I_{info}-I_{control}')
    title('Diferencia por acoplamiento')
    % Figuras 03 y 04: paneles completos de ambos escenarios.
    if c==1
        archivoFigura='fig03_escenario_desinformacion.png';
    else
        archivoFigura='fig04_escenario_escepticismo.png';
    end
    exportgraphics(gcf,fullfile(carpetaFiguras,archivoFigura),'Resolution',300);
    fprintf('Figura exportada: %s\n',archivoFigura);
end

%% 3. Exportar resultados numericos
if ~isempty(filas)
    tabla=cell2table(filas,'VariableNames',{'Escenario','Tratamiento','R00','R01','R02', ...
        'I_Final','I_Prom_Ultimo10','I_Pico','Error_RK4_vs_ode15s', ...
        'Error_Conservacion','Minimo'});
    disp(tabla)
    writetable(tabla,'MN_Grupo8_resultados_escenarios.csv');
    fprintf('\nGuardado: MN_Grupo8_resultados_escenarios.csv\n');
end
fprintf('NOTA: Las comparaciones son exploratorias, no equivalen a demostrar extincion o endemia asintotica.\n');

%% FUNCIONES LOCALES
function [t,X,meta]=rk4_integrar(x0,p,T,h)
    n=ceil(T/h); h=T/n;
    t=linspace(0,T,n+1);
    X=zeros(6,n+1); X(:,1)=x0;
    meta.valido=true;
    for k=1:n
        a=f_umzsirs(X(:,k),p);
        b=f_umzsirs(X(:,k)+h*a/2,p);
        c=f_umzsirs(X(:,k)+h*b/2,p);
        d=f_umzsirs(X(:,k)+h*c,p);
        X(:,k+1)=X(:,k)+h*(a+2*b+2*c+d)/6;
        y=X(:,k+1);
        if any(~isfinite(y)) || min(y)<-1e-10 || 1-p.K*y(5)<=0 || 1+y(3)<=0
            meta.valido=false;
            X=X(:,1:k+1); t=t(1:k+1);
            break
        end
    end
    meta.cons=max([max(abs(sum(X(1:3,:),1)-1)), max(abs(sum(X(4:6,:),1)-1))]);
    meta.minimo=min(X,[],'all');
end

function dx=f_umzsirs(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5); R=x(6);
    a=p.b2/(1-p.K*I); b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end

function s=si_no(cond)
    if cond, s='SI'; else, s='NO'; end
end
