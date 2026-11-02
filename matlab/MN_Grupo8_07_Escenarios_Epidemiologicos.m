%% Grupo 8 - UMZ-SIRS | 07. Escenarios epidemiologicos (ENTREGA FINAL)
% Metodo numerico principal: RK4 implementado por el grupo (h fijo).
% Validacion independiente: ode15s con tolerancia ESTRICTA (no se presenta como
% metodo propio) y ode15s con tolerancia laxa para medir la diferencia entre
% referencias. Ademas se compara con el equilibrio endemico ANALITICO.
% Escenarios I y II: parametros consignados en el articulo base (Bulai et al., 2024).
% Condiciones iniciales y tau_final: decisiones del grupo, no del articulo.
% Salidas: ../resultados/MN_Grupo8_resultados_escenarios.csv
%          ../figuras/fig03_escenario_desinformacion.png y fig04_escenario_escepticismo.png
% Nota: guardar en UTF-8 (los textos de las figuras llevan tildes).
clear; clc; close all; format long g

carpetaFiguras='../figuras'; carpetaRes='../resultados';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Escenarios y condiciones iniciales
p.mu1=5.2143; p.mu2=5.2143; p.gamma=0.1667;
p.eta=0.1; p.K=0.9; p.eps=1/50;
tauFinal=30; T=tauFinal/p.eps;
h=0.1;                     % paso fijo en tiempo rapido
I0=0.02; S0=1-I0; R0=0;

casos(1).nombre='I - Desinformación'; casos(1).tag='escenario I (desinformación)';
casos(1).b1=6.7041; casos(1).b2=5.9041; casos(1).beta=4.8429;
casos(1).informacion=[0.8;0.2;0];
casos(2).nombre='II - Escepticismo'; casos(2).tag='escenario II (escepticismo)';
casos(2).b1=5.9041; casos(2).b2=6.7041; casos(2).beta=5.919;
casos(2).informacion=[0.8;0;0.2];
infoControl=[1;0;0];       % control: U=1, M=Z=0 es un estado invariante exacto

optEstricta=odeset('RelTol',1e-12,'AbsTol',1e-14,'MaxStep',1);
optLaxa=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',1);

%% 2. Integracion y comparacion con control
filas={};
for c=1:numel(casos)
    q=p; q.b1=casos(c).b1; q.b2=casos(c).b2; q.beta=casos(c).beta;
    % Umbrales tomados del borrador del grupo; verificar contra Bulai et al. (2024).
    R00=q.beta/(q.gamma+q.mu2);
    R01=R00*(2*q.b1-q.mu1)/q.b1;
    R02=R00*q.b2/(2*q.b2-q.mu1);
    eq=equilibrio_analitico(q);
    Iend=@(be) max(0,(1-(q.gamma+q.mu2)./be)./(1+q.gamma/(q.eta+q.mu2)));
    IeqControl=Iend(q.beta);
    if c==1, IeqAcopl=eq.M.I; okRama=eq.M.ok; else, IeqAcopl=eq.Z.I; okRama=eq.Z.ok; end
    fprintf('\n===== %s =====\n',casos(c).nombre);
    fprintf('R00=%.6f, R01=%.6f, R02=%.6f | I* analitico: control=%.6g, acoplado=%.6g (rama coherente: %d)\n', ...
        R00,R01,R02,IeqControl,IeqAcopl,okRama);

    x0control=[infoControl;S0;I0;R0];
    x0acoplado=[casos(c).informacion;S0;I0;R0];
    [tC,XC,mC]=rk4_integrar(x0control,q,T,h);
    [tA,XA,mA]=rk4_integrar(x0acoplado,q,T,h);
    if ~mC.valido || ~mA.valido
        warning('Caso %d: RK4 dejo dominio admisible. Se omite conclusion fisica.',c);
        continue
    end
    % Dos referencias ode15s: estricta (se usa como referencia) y laxa (para medir
    % cuanto cambia la referencia al aflojar la tolerancia).
    solCe=ode15s(@(t,x)f_umzsirs(x,q),[0 T],x0control,optEstricta);
    solAe=ode15s(@(t,x)f_umzsirs(x,q),[0 T],x0acoplado,optEstricta);
    solCl=ode15s(@(t,x)f_umzsirs(x,q),[0 T],x0control,optLaxa);
    solAl=ode15s(@(t,x)f_umzsirs(x,q),[0 T],x0acoplado,optLaxa);
    errC=max(abs(XC-deval(solCe,tC)),[],'all');
    errA=max(abs(XA-deval(solAe,tA)),[],'all');
    difC=max(abs(deval(solCe,tC)-deval(solCl,tC)),[],'all');
    difA=max(abs(deval(solAe,tA)-deval(solAl,tA)),[],'all');
    % Promedio en el ultimo 10 % del tiempo lento: no confundir con limite t->infinito.
    ultC=tC>=0.9*T; ultA=tA>=0.9*T;
    IfinC=XC(5,end); IfinA=XA(5,end);
    IpromC=mean(XC(5,ultC)); IpromA=mean(XA(5,ultA));
    picoC=max(XC(5,:)); picoA=max(XA(5,:));
    fprintf('RK4, h=%g, tau final=%g: %d pasos por trayectoria\n',h,tauFinal,numel(tC)-1);
    fprintf('Control: I_final=%.8f (analitico %.8f), I_pico=%.8f\n',IfinC,IeqControl,picoC);
    fprintf('Informacion: I_final=%.8f (analitico %.8f), I_pico=%.8f\n',IfinA,IeqAcopl,picoA);
    fprintf('Error max RK4 vs ode15s estricta: control %.3e | informacion %.3e\n',errC,errA);
    fprintf('Diferencia entre referencias (estricta vs laxa): control %.3e | informacion %.3e\n',difC,difA);
    fprintf('Conservacion max: control %.3e | informacion %.3e\n',mC.cons,mA.cons);

    filas(end+1,:)={c,'Control',R00,R01,R02,IfinC,IpromC,picoC,IeqControl, ...
        IfinC-IeqControl,errC,difC,mC.cons,mC.minimo}; %#ok<AGROW>
    filas(end+1,:)={c,'Con informacion',R00,R01,R02,IfinA,IpromA,picoA,IeqAcopl, ...
        IfinA-IeqAcopl,errA,difA,mA.cons,mA.minimo}; %#ok<AGROW>

    %% Figura 2x2: (a) I con/sin informacion, (b) diferencia, (c) 6 variables con
    % informacion, (d) 6 variables del control. Paneles (c) y (d) en tau logaritmico
    % para que la capa rapida (tau pequeno) y la lenta se vean a la vez.
    fig=figure('Color','w','Name',casos(c).nombre,'Position',[40 40 1250 800]);
    nom={'U','M','Z','S','I','R'}; col=lines(6);
    subplot(2,2,1)
    plot(q.eps*tC,XC(5,:),'--','LineWidth',1.5); hold on
    plot(q.eps*tA,XA(5,:),'LineWidth',1.5); grid on
    xlabel('Tiempo lento \tau'); ylabel('Fracción infectada I')
    title('(a) Infección con y sin difusión informativa')
    legend('Control (sin información)','Con información','Location','best')
    subplot(2,2,2)
    plot(q.eps*tA,XA(5,:)-XC(5,:),'LineWidth',1.5); yline(0,'--'); grid on
    xlabel('Tiempo lento \tau'); ylabel('I_{info} - I_{control}')
    title('(b) Efecto del acoplamiento informativo')
    subplot(2,2,3); hold on
    for k=1:6
        if k<=3, est='--'; else, est='-'; end
        semilogx(q.eps*tA(2:end),XA(k,2:end),est,'Color',col(k,:),'LineWidth',1.4);
    end
    set(gca,'XScale','log'); grid on; xlim([1e-3 tauFinal]); ylim([-0.02 1.02]);
    xlabel('Tiempo lento \tau (escala log.)'); ylabel('Proporción')
    title(['(c) Con información: ' casos(c).tag])
    legend(nom,'NumColumns',3,'Location','best')
    subplot(2,2,4); hold on
    for k=1:6
        if k<=3, est='--'; else, est='-'; end
        semilogx(q.eps*tC(2:end),XC(k,2:end),est,'Color',col(k,:),'LineWidth',1.4);
    end
    set(gca,'XScale','log'); grid on; xlim([1e-3 tauFinal]); ylim([-0.02 1.02]);
    xlabel('Tiempo lento \tau (escala log.)'); ylabel('Proporción')
    title('(d) Control: sin difusión informativa (U=1, M=Z=0)')
    legend(nom,'NumColumns',3,'Location','best')
    sgtitle({sprintf('%s: (b_1,b_2,\\beta)=(%g, %g, %g), \\epsilon=%g', ...
        casos(c).nombre,q.b1,q.b2,q.beta,q.eps), ...
        sprintf('(U,M,Z)_0=(%g, %g, %g); líneas discontinuas: capa informativa', ...
        casos(c).informacion(1),casos(c).informacion(2),casos(c).informacion(3))},'FontSize',12);
    if c==1, archivoFigura='fig03_escenario_desinformacion.png';
    else, archivoFigura='fig04_escenario_escepticismo.png'; end
    exportgraphics(fig,fullfile(carpetaFiguras,archivoFigura),'Resolution',300);
    fprintf('Figura exportada: %s\n',archivoFigura);
end

%% 3. Exportar resultados numericos
if ~isempty(filas)
    tabla=cell2table(filas,'VariableNames',{'Escenario','Tratamiento','R00','R01','R02', ...
        'I_Final','I_Prom_Ultimo10','I_Pico','I_Equilibrio_Analitico','Dif_Final_vs_Analitico', ...
        'Error_RK4_vs_ode15s','Dif_entre_referencias','Error_Conservacion','Minimo'});
    disp(tabla)
    writetable(tabla,fullfile(carpetaRes,'MN_Grupo8_resultados_escenarios.csv'));
    fprintf('\nGuardado: MN_Grupo8_resultados_escenarios.csv\n');
end
fprintf('NOTA: Las comparaciones son de horizonte finito; no demuestran extincion o endemia asintotica.\n');

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
    if 1-p.K*I<=0 || 1+Z<=0, dx=NaN(6,1); return, end
    a=p.b2/(1-p.K*I); b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end

function eq=equilibrio_analitico(p)
    % Equilibrio de la capa lenta: S*=(gamma+mu2)/beta_ef, R*=gamma I*/(eta+mu2),
    % S*+I*+R*=1  =>  I*=(1-(gamma+mu2)/beta_ef)/(1+gamma/(eta+mu2)), I*>=0.
    Iend=@(be) max(0,(1-(p.gamma+p.mu2)./be)./(1+p.gamma/(p.eta+p.mu2)));
    % Rama M (Z=0): U=mu1/b1, M=1-mu1/b1 (si b1>mu1), beta_ef=beta(1+M).
    Ms=max(0,1-p.mu1/p.b1); I_M=Iend(p.beta*(1+Ms));
    eq.M=struct('I',I_M,'U',1-Ms,'M',Ms,'Z',0, ...
        'ok',(p.b1>p.mu1)&&(p.b2/(1-p.K*I_M)<=p.b1));   % Z no puede invadir
    % Rama Z (M=0): U=mu1(1-K I)/b2, Z=1-U, beta_ef=beta/(1+Z); I* unica (biseccion).
    Zf=@(I) max(0,1-p.mu1*(1-p.K*I)/p.b2);
    g=@(I) I-Iend(p.beta/(1+Zf(I)));
    lo=0; hi=1;
    for k=1:60
        mid=(lo+hi)/2;
        if g(mid)<0, lo=mid; else, hi=mid; end
    end
    I_Z=(lo+hi)/2; Zz=Zf(I_Z);
    eq.Z=struct('I',I_Z,'U',1-Zz,'M',0,'Z',Zz, ...
        'ok',(p.b2/(1-p.K*I_Z)>p.mu1)&&(p.b1*(1-p.K*I_Z)/p.b2<=1));   % M no puede invadir
    eq.I0=Iend(p.beta);
    if eq.M.ok, eq.regimen='M'; eq.I=eq.M.I;
    elseif eq.Z.ok, eq.regimen='Z'; eq.I=eq.Z.I;
    elseif p.b1<=p.mu1 && p.b2/(1-p.K*eq.I0)<=p.mu1, eq.regimen='ninguno'; eq.I=eq.I0;
    else, eq.regimen='indefinido'; eq.I=NaN;
    end
end
