%% Grupo 8 - UMZ-SIRS | 11. Costo computacional frente a epsilon (NUEVO)
% Pregunta: cuanto cuesta alcanzar una precision dada cuando 1/epsilon crece?
% Metodos: RK4 (paso fijo), Euler atras (solucion directa del sistema lineal,
% paso fijo), ode45 (explicito adaptativo) y ode15s (implicito adaptativo).
% Dos conjuntos: base exploratorio (mu=1) y escenario I de Bulai et al. (mu=5.2).
% Referencia: ode15s con tolerancia estricta. Para cada (conjunto, epsilon, metodo)
% se prueban varios parametros (h o RelTol) y la figura usa el de MENOR costo que
% alcanza el error objetivo.
% AJUSTE: si tarda demasiado, poner modoRapido=true (epsilon hasta 1e-4).
% Salidas: ../resultados/MN_Grupo8_resultados_costo_epsilon.csv
%          ../figuras/fig07_costo_vs_epsilon.png
% Nota: guardar en UTF-8 (los textos de las figuras llevan tildes).
clear; clc; close all; format short g

carpetaFiguras='../figuras'; carpetaRes='../resultados';

% LETRA LEGIBLE AL ANCHO DE COLUMNA DEL PAPER

set(groot,'defaultAxesFontSize',13,'defaultLegendFontSize',11,'defaultTextFontSize',13);

if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Configuracion
modoRapido=false;
if modoRapido, epsilons=10.^(-1:-1:-4); else, epsilons=10.^(-1:-1:-6); end
errObjetivo=1e-3;
maxPasosExpl=1e6; maxPasosImpl=1e5;
rtolsODE=[1e-3 1e-5 1e-7];
optRef=odeset('RelTol',1e-11,'AbsTol',1e-13);

opt.tolResidual=1e-10; opt.tolGS=1e-6; opt.tolLin=1e-11;
opt.maxNewton=50; opt.maxGS=50; opt.tolPos=1e-10; opt.tolCons=1e-8;
opt.minH=1e-8; opt.maxRechazos=50; opt.reducir=false; opt.solver='directo';

c1.nombre='Base exploratorio (\mu=1)';
c1.p=struct('b1',1.1,'b2',0.92,'mu1',1,'K',0.9,'beta',6,'gamma',1.2,'eta',0.7,'mu2',1,'eps',0);
c1.x0=[0.7;0.2;0.1;0.8;0.2;0]; c1.tauFinal=1;
c1.hRK4=[0.5 1 2]; c1.hEA=[0.5 2 10 50];
c2.nombre='Escenario I de Bulai (\mu=5,2)';
c2.p=struct('b1',6.7041,'b2',5.9041,'mu1',5.2143,'K',0.9,'beta',4.8429,'gamma',0.1667, ...
    'eta',0.1,'mu2',5.2143,'eps',0);
c2.x0=[0.8;0.2;0;0.98;0.02;0]; c2.tauFinal=5;
c2.hRK4=[0.1 0.2 0.4]; c2.hEA=[0.1 1 10 50];
conjuntos={c1,c2};

%% 2. Barrido
filas={};
for ic=1:2
    cj=conjuntos{ic};
    for ie=1:numel(epsilons)
        p=cj.p; p.eps=epsilons(ie); T=cj.tauFinal/p.eps; x0=cj.x0;
        fun=@(t,x) modelo(x,p);
        fprintf('\n== %s | eps=%g | T=%g ==\n',strrep(cj.nombre,'\mu','mu'),p.eps,T);
        ref=ode15s(fun,[0 T],x0,optRef);
        % Malla uniforme + logaritmica: incluye el transitorio rapido aunque T sea enorme.
        tg=unique(min([linspace(0,T,2001), logspace(-2,log10(T),2001)],T)); Yref=deval(ref,tg);

        % RK4 de paso fijo
        for h=cj.hRK4
            filas(end+1,:)=correr_fijo(ic,cj.nombre,p.eps,T,'RK4',h,x0,p,opt,ref,maxPasosExpl); %#ok<AGROW>
        end
        % Euler atras (solucion directa del sistema lineal)
        for h=cj.hEA
            filas(end+1,:)=correr_fijo(ic,cj.nombre,p.eps,T,'Euler atras',h,x0,p,opt,ref,maxPasosImpl); %#ok<AGROW>
        end
        % ode45 y ode15s con varias tolerancias
        for rt=rtolsODE
            for im=1:2
                if im==1, nom='ode45'; sv=@ode45; else, nom='ode15s'; sv=@ode15s; end
                try
                    reloj=tic;
                    sol=sv(fun,[0 T],x0,odeset('RelTol',rt,'AbsTol',1e-2*rt));
                    seg=toc(reloj);
                    err=max(abs(deval(sol,tg)-Yref),[],'all');
                    st=sol.stats; npd=0;
                    if isfield(st,'npds'), npd=st.npds; end
                    fila={ic,cj.nombre,p.eps,T,nom,rt,'OK',err,st.nsteps,st.nfailed, ...
                        st.nfevals,npd,seg};
                catch ME
                    fila={ic,cj.nombre,p.eps,T,nom,rt,['FALLO: ' ME.message],NaN,NaN,NaN,NaN,NaN,NaN};
                end
                filas(end+1,:)=fila; %#ok<AGROW>
            end
        end
    end
end
tabla=cell2table(filas,'VariableNames',{'Conjunto','NombreConjunto','Epsilon','T', ...
    'Metodo','Parametro','Estado','ErrorMaxRef','Pasos','Rechazados','EvalF', ...
    'Jacobianas','Segundos'});
writetable(tabla,fullfile(carpetaRes,'MN_Grupo8_resultados_costo_epsilon.csv'));
fprintf('\nGuardado: MN_Grupo8_resultados_costo_epsilon.csv\n');

%% 3. Figura: costo minimo (evaluaciones de f) que alcanza el error objetivo
metodos={'RK4','Euler atras','ode45','ode15s'};
etiq={'RK4','Euler atrás (directo)','DP45','NDF'};
col=lines(4); mk={'o-','s-','^-','d-'};
fig=figure('Color','w','Name','Costo vs epsilon','Position',[50 50 1200 480]);
for ic=1:2
    subplot(1,2,ic); hold on;
    hl=gobjects(1,4);
    for im=1:4
        costo=NaN(1,numel(epsilons));
        for ie=1:numel(epsilons)
            m=tabla.Conjunto==ic & tabla.Epsilon==epsilons(ie) & strcmp(tabla.Metodo,metodos{im}) ...
                & strcmp(tabla.Estado,'OK') & tabla.ErrorMaxRef<=errObjetivo;
            if any(m), costo(ie)=min(tabla.EvalF(m)); end
        end
        hl(im)=plot(1./epsilons,costo,mk{im},'Color',col(im,:),'LineWidth',1.4,'MarkerSize',6);
    end
    set(gca,'XScale','log','YScale','log'); grid on;
    xlabel('1/\epsilon'); ylabel('Evaluaciones de f (mínimo para error \leq 10^{-3})');
    title(conjuntos{ic}.nombre); legend(hl,etiq,'Location','northwest');
end
exportgraphics(fig,fullfile(carpetaFiguras,'fig07_costo_vs_epsilon.png'),'Resolution',300);
fprintf('Figura exportada: fig07_costo_vs_epsilon.png\n');
fprintf('Las curvas sin punto indican que ningun parametro probado alcanzo el error objetivo.\n');
fprintf('EvalF no incluye Jacobianas ni resolucion de sistemas lineales (ver columna Jacobianas).\n');

%% FUNCIONES LOCALES
function fila=correr_fijo(ic,nombre,eps0,T,metodo,h,x0,p,opt,ref,maxPasos)
    if T/h>maxPasos
        fila={ic,nombre,eps0,T,metodo,h,'OMITIDO (demasiados pasos)',NaN,NaN,NaN,NaN,NaN,NaN};
        return
    end
    reloj=tic;
    if strcmp(metodo,'RK4'), [t,X,info]=explicito(x0,p,T,h,4,opt);
    else, [t,X,info]=euler_atras(x0,p,T,h,opt); end
    seg=toc(reloj);
    if strcmp(info.estado,'OK')
        err=max(abs(X-deval(ref,min(max(t,0),T))),[],'all');
    else
        err=NaN;
    end
    if ~isfield(info,'jac'), info.jac=0; end
    fila={ic,nombre,eps0,T,metodo,h,info.estado,err,info.aceptados,info.rechazados, ...
        info.evalf,info.jac,seg};
end

function [t,X,info]=explicito(x0,p,T,hBase,orden,opt)
    cap=1024; tb=zeros(1,cap); Xb=zeros(6,cap); n=1; Xb(:,1)=x0;
    h=hBase; seguidos=0;
    info=struct('estado','OK','aceptados',0,'rechazados',0,'evalf',0);
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
            info.aceptados=info.aceptados+1;
            h=min(hBase,1.5*h); seguidos=0;
        else
            info.rechazados=info.rechazados+1; seguidos=seguidos+1;
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

function [t,X,info]=euler_atras(x0,p,T,hBase,opt)
    cap=1024; tb=zeros(1,cap); Xb=zeros(6,cap); n=1; Xb(:,1)=x0;
    h=hBase; seguidos=0;
    info=struct('estado','OK','aceptados',0,'rechazados',0,'newton',0,'gs',0, ...
        'jac',0,'evalf',0,'rhoMax',0);
    while T-tb(n)>1e-10*max(1,T)
        h=min(h,T-tb(n));
        [xn,ok,causa,s]=paso_implicito(Xb(:,n),p,h,opt);
        info.newton=info.newton+s.newton; info.gs=info.gs+s.gs;
        info.jac=info.jac+s.jac; info.evalf=info.evalf+s.evalf;
        if ok && ~(all(xn>=-opt.tolPos) && abs(sum(xn(1:3))-1)<opt.tolCons ...
                && abs(sum(xn(4:6))-1)<opt.tolCons)
            ok=false; causa='positividad';
        end
        if ok
            n=n+1; if n>cap, cap=2*cap; tb(cap)=0; Xb(6,cap)=0; end
            tb(n)=tb(n-1)+h; Xb(:,n)=xn;
            info.aceptados=info.aceptados+1;
            h=min(hBase,1.5*h); seguidos=0;
        else
            info.rechazados=info.rechazados+1; seguidos=seguidos+1;
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
