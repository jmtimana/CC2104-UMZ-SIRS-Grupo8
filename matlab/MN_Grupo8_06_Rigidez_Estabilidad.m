%% GRUPO 8 - Experimento 06: rigidez y estabilidad UMZ-SIRS
% Valores de epsilon y pasos temporales. Comparacion con ode15s.
% Ejecutar el archivo completo en MATLAB; guardar luego como .mlx si procede.
clear; clc; close all; format short g

%% Exportacion automatica para el paper
carpetaFiguras = 'figuras';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end

%% 1. Modelo, parametros y horizonte epidemiologico comun
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1;
x0=[0.7;0.2;0.1;0.8;0.2;0];
epsilons=[1e-1 1e-2 1e-3];
pasos=[0.5 1 2 4];
tauFinal=1;
tolPos=1e-10; tolCons=1e-8;
opt.tolNewton=1e-8; opt.tolResidual=1e-10;
opt.tolGS=1e-6; opt.tolResiduoLineal=1e-11;
opt.maxNewton=50; opt.maxGS=50;
opt.tolPos=tolPos; opt.tolCons=tolCons;
opt.minH=1e-8; opt.maxRechazos=50;

%% 2. Simulaciones con los mismos horizontes lentos
Resultados=table();
figure('Color','w','Name','Rigidez: error vs paso');
for ie=1:numel(epsilons)
    p.eps=epsilons(ie); T=tauFinal/p.eps;
    fun=@(t,x) modelo(x,p);
    opciones=odeset('RelTol',1e-11,'AbsTol',1e-13);
    referencia=ode15s(fun,[0 T],x0,opciones);
    J=jacobiana_f(x0,p);
    aut=eig(J);
    reales=abs(real(aut)); reales=reales(reales>1e-10);
    ratio=max(reales)/min(reales);
    fprintf('\n==== epsilon=%g, T=%g, razon espectral inicial=%g ====\n',p.eps,T,ratio);

    for ih=1:numel(pasos)
        h=pasos(ih);
        for im=1:3
            metodos={'Euler explicito','RK4','Euler atras'};
            nombre=metodos{im};
            tic;
            try
                if im<3
                    [t,X,nf,estado]=integrar_explicito(x0,p,T,h,im,tolPos);
                    rechazados=0; nNewton=0; nGS=0;
                else
                    [t,X,info]=integrar_implicito(x0,p,T,h,opt);
                    nf=info.evalf; rechazados=info.rechazados;
                    nNewton=info.newton; nGS=info.gs;
                    estado='OK';
                end
                segundos=toc;
                if strcmp(estado,'OK')
                    Xref=deval(referencia,t);
                    errorRef=max(abs(X-Xref),[],'all');
                    cons=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
                    minimo=min(X,[],'all');
                    if cons>tolCons || minimo < -tolPos
                        estado='VIOLA RESTRICCIONES';
                    end
                else
                    errorRef=NaN; cons=NaN; minimo=NaN;
                end
                pasosHechos=numel(t)-1;
            catch ME
                segundos=toc;
                estado=['FALLO: ' ME.message];
                errorRef=NaN; cons=NaN; minimo=NaN;
                nf=NaN; rechazados=NaN; nNewton=NaN; nGS=NaN;
                pasosHechos=NaN;
            end
            fila=table(p.eps,T,h,{nombre},errorRef,cons,minimo,pasosHechos, ...
                rechazados,nf,nNewton,nGS,segundos,{estado},ratio, ...
                'VariableNames',{'Epsilon','T','h','Metodo','ErrorMaxRef', ...
                'ErrorConservacion','Minimo','Pasos','Rechazados','EvalF', ...
                'IterNewton','IterGS','Segundos','Estado','RazonEspectral'});
            Resultados=[Resultados;fila]; %#ok<AGROW>
        end
    end
    subplot(1,3,ie); hold on
    for im=1:3
        mascara=Resultados.Epsilon==p.eps & strcmp(Resultados.Metodo,metodos{im}) ...
            & strcmp(Resultados.Estado,'OK') & Resultados.ErrorMaxRef>0;
        loglog(Resultados.h(mascara),Resultados.ErrorMaxRef(mascara),'-o','LineWidth',1.1);
    end
    xlabel('Paso nominal h'); ylabel('Error maximo vs ode15s');
    title(sprintf('epsilon = %g',p.eps)); grid on;
    legend(metodos,'Location','best');
end

% Figura 02: tres valores de epsilon en una unica imagen.
exportgraphics(gcf,fullfile(carpetaFiguras, ...
    'fig02_rigidez_estabilidad.png'),'Resolution',300);
fprintf('Figura exportada: fig02_rigidez_estabilidad.png\n');

%% 3. Guardar resultados para tablas del paper
disp(Resultados);
writetable(Resultados,'MN_Grupo8_resultados_rigidez.csv');
fprintf('\nGuardado: MN_Grupo8_resultados_rigidez.csv\n');
fprintf('Nota: el paso implicito puede reducir h por rechazo. Comparar Pasos, no solo h.\n');
fprintf('La razon espectral es un indicador local, no prueba por si sola rigidez.\n');

%% FUNCIONES LOCALES
function [t,X,nf,estado]=integrar_explicito(x0,p,T,h,metodo,tolPos)
    n=ceil(T/h); t=zeros(1,n+1); X=zeros(6,n+1); X(:,1)=x0;
    nf=0; estado='OK';
    for k=1:n
        hk=min(h,T-t(k)); x=X(:,k);
        if metodo==1
            nuevo=x+hk*modelo(x,p); nf=nf+1;
        else
            k1=modelo(x,p); k2=modelo(x+hk*k1/2,p);
            k3=modelo(x+hk*k2/2,p); k4=modelo(x+hk*k3,p);
            nuevo=x+hk*(k1+2*k2+2*k3+k4)/6; nf=nf+4;
        end
        t(k+1)=t(k)+hk; X(:,k+1)=nuevo;
        if any(~isfinite(nuevo)) || any(abs(nuevo)>1e8) || ...
                nuevo(5)*p.K>=1 || nuevo(3)<=-1 || any(nuevo<-tolPos)
            t=t(1:k+1); X=X(:,1:k+1);
            estado='INADMISIBLE/INESTABLE'; return
        end
    end
end

function [t,X,info]=integrar_implicito(x0,p,T,hBase,opt)
    t=0; X=x0; h=hBase;
    info.rechazados=0;info.newton=0;info.gs=0;info.evalf=0;
    rSeguidos=0;
    while T-t(end)>1e-10*max(1,T)
        h=min(h,T-t(end)); anterior=X(:,end);
        [nuevo,ok,nn,ng,nf]=paso_implicito(anterior,p,h,opt);
        info.newton=info.newton+nn; info.gs=info.gs+ng;info.evalf=info.evalf+nf;
        if ok && all(nuevo>=-opt.tolPos) && ...
             abs(sum(nuevo(1:3))-1)<opt.tolCons && ...
             abs(sum(nuevo(4:6))-1)<opt.tolCons
            t(end+1)=t(end)+h; X(:,end+1)=nuevo; %#ok<AGROW>
            h=min(hBase,1.5*h);rSeguidos=0;
        else
            info.rechazados=info.rechazados+1; rSeguidos=rSeguidos+1;
            h=h/2;
            if h<opt.minH || rSeguidos>opt.maxRechazos
                error('No converge Newton/GS incluso tras reducir el paso');
            end
        end
    end
end

function [x,ok,nn,ng,nf]=paso_implicito(x0,p,h,opt)
    x=x0;ok=false;ng=0;nf=0;nn=0;
    for k=1:opt.maxNewton
        nn=k; fx=modelo(x,p);nf=nf+1;
        G=x-x0-h*fx;
        if norm(G,inf)<opt.tolResidual
            ok=true;return
        end
        W=eye(6)-h*jacobiana_f(x,p);
        if any(~isfinite(W),'all') || any(abs(diag(W))<1e-13),return,end
        B=-(tril(W)\triu(W,1));
        if max(abs(eig(B)))>=1,return,end
        [delta,bien,iter]=gs(W,-G,opt);
        ng=ng+iter;
        if ~bien,return,end
        alpha=1;pasoValido=false;
        while alpha>=1/1024
            candidato=x+alpha*delta;
            if all(isfinite(candidato)) && all(candidato>=-opt.tolPos) ...
                  && 1-p.K*candidato(5)>0 && 1+candidato(3)>0
                v=candidato-x0-h*modelo(candidato,p);nf=nf+1;
                if norm(v,inf)<norm(G,inf)
                    x=candidato;pasoValido=true;break
                end
            end
            alpha=alpha/2;
        end
        if ~pasoValido,return,end
    end
end

function [y,ok,k]=gs(A,b,opt)
    n=numel(b);y=zeros(n,1);ok=false;
    for k=1:opt.maxGS
        ant=y;
        for i=1:n
            y(i)=(b(i)-A(i,1:i-1)*y(1:i-1)-A(i,i+1:end)*ant(i+1:end))/A(i,i);
        end
        if norm(y-ant,inf)<opt.tolGS && norm(A*y-b,inf)<opt.tolResiduoLineal
            ok=true;return
        end
    end
end

function dx=modelo(x,p)
    U=x(1);M=x(2);Z=x(3);S=x(4);I=x(5);R=x(6);
    if 1-p.K*I<=0 || 1+Z<=0,error('Denominador no positivo'),end
    a=p.b2/(1-p.K*I);b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end

function J=jacobiana_f(x,p)
    U=x(1);M=x(2);Z=x(3);S=x(4);I=x(5);
    a=p.b2/(1-p.K*I); da=p.b2*p.K/(1-p.K*I)^2;
    b=p.beta*(1+M)/(1+Z); bm=p.beta/(1+Z);
    bz=-p.beta*(1+M)/(1+Z)^2;e=p.eps;
    J=[-p.b1*M-a*Z-p.mu1,-p.b1*U,-a*U,0,-da*U*Z,0;
        p.b1*M,p.b1*U-p.mu1,0,0,0,0;
        a*Z,0,a*U-p.mu1,0,da*U*Z,0;
        0,-e*bm*S*I,-e*bz*S*I,e*(-b*I-p.mu2),-e*b*S,e*p.eta;
        0,e*bm*S*I,e*bz*S*I,e*b*I,e*(b*S-p.gamma-p.mu2),0;
        0,0,0,0,e*p.gamma,-e*(p.eta+p.mu2)];
end
