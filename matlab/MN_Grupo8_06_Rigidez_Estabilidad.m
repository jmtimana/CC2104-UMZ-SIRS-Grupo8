%% GRUPO 8 - Experimento 06: rigidez y estabilidad UMZ-SIRS (ENTREGA FINAL)
% Ejecutar Run con matlab/ como directorio de trabajo.
% Compara, para varios epsilon y pasos nominales h:
%   - Euler explicito y RK4: politica NOMINAL (sin reducir h) y con REDUCCION de paso
%   - RKF45 adaptativo con tope h_max = h nominal
%   - Euler atras con Gauss-Seidel y con solucion directa (control), ambas politicas
% Ademas: h maximo estable por biseccion (horizonte largo comun, para no mezclar
% rigidez con longitud del horizonte), rho(T_GS) frente a h y trayectoria del paso.
% Referencia: ode15s con tolerancias estrictas (se valida en el archivo 09).
% Salidas: ../resultados/MN_Grupo8_resultados_rigidez.csv
%          ../resultados/MN_Grupo8_resultados_hmax_estable.csv
%          ../resultados/MN_Grupo8_resultados_rho_gs.csv
%          ../figuras/fig02_rigidez_estabilidad.png, fig08_rho_gauss_seidel.png,
%          fig09_paso_efectivo.png
% Nota: guardar en UTF-8 (los textos de las figuras llevan tildes).
clear; clc; close all; format short g

carpetaFiguras='../figuras'; carpetaRes='../resultados';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Modelo, parametros y horizonte epidemiologico comun
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1;
x0=[0.7;0.2;0.1;0.8;0.2;0];
epsilons=[1e-1 1e-2 1e-3];
pasos=[0.5 1 2 4];
tauFinal=1;
Tbis=1000;                          % horizonte (tiempo rapido) para el h maximo estable
opt.tolResidual=1e-10; opt.tolGS=1e-6; opt.tolLin=1e-11;
opt.maxNewton=50; opt.maxGS=50; opt.tolPos=1e-10; opt.tolCons=1e-8;
opt.minH=1e-8; opt.maxRechazos=50; opt.hmin=1e-10; opt.maxIntentos=500000;
opt.rtol=1e-8; opt.atol=1e-10;      % tolerancias de RKF45
opt.reducir=false; opt.solver='gs';
optRef=odeset('RelTol',1e-11,'AbsTol',1e-13);

% Configuraciones: metodo, politica, solver lineal
cfg={'Euler explicito','nominal','-';
     'Euler explicito','reduccion','-';
     'RK4','nominal','-';
     'RK4','reduccion','-';
     'RKF45','adaptativo','-';
     'Euler atras','nominal','GS';
     'Euler atras','reduccion','GS';
     'Euler atras','nominal','directo';
     'Euler atras','reduccion','directo'};

%% 2. Simulaciones con los mismos horizontes lentos
filas={}; filasH={}; trayectorias={};
for ie=1:numel(epsilons)
    p.eps=epsilons(ie); T=tauFinal/p.eps;
    fun=@(t,x) modelo(x,p);
    ref=ode15s(fun,[0 T],x0,optRef);
    J=jacobiana_f(x0,p);
    reales=abs(real(eig(J))); reales=reales(reales>1e-10);
    ratio=max(reales)/min(reales);
    % Autovalor dominante a lo largo de la referencia (para comparar con 2/|lambda|)
    tm=linspace(0,Tbis,200); Xm=deval(ode15s(fun,[0 Tbis],x0,optRef),tm);
    lmaxTray=0; for k=1:size(Xm,2), lmaxTray=max(lmaxTray,max(abs(real(eig(jacobiana_f(Xm(:,k),p)))))); end
    lmaxFinal=max(abs(real(eig(jacobiana_f(Xm(:,end),p)))));
    fprintf('\n==== epsilon=%g, T=%g, razon espectral inicial=%g, |lambda|max final=%.4g ====\n', ...
        p.eps,T,ratio,lmaxFinal);

    for ih=1:numel(pasos)
        h=pasos(ih);
        for ic=1:size(cfg,1)
            [fila,t]=ejecutar(cfg{ic,1},cfg{ic,2},cfg{ic,3},x0,p,T,h,opt,ref,ratio);
            filas(end+1,:)=fila; %#ok<AGROW>
            if p.eps==1e-2 && h==4 && ~strcmp(cfg{ic,2},'nominal')
                trayectorias{end+1}={cfg{ic,1},cfg{ic,3},t}; %#ok<AGROW>
            end
        end
    end

    % h maximo estable (politica nominal, horizonte largo comun)
    metH={'Euler explicito','-'; 'RK4','-'; 'Euler atras','GS'; 'Euler atras','directo'};
    for im=1:size(metH,1)
        hm=biseccion_hmax(metH{im,1},metH{im,2},x0,p,Tbis,opt);
        if strcmp(metH{im,1},'Euler explicito'), lim=2/lmaxFinal;
        elseif strcmp(metH{im,1},'RK4'), lim=2.78/lmaxFinal; else, lim=NaN; end
        filasH(end+1,:)={p.eps,metH{im,1},metH{im,2},hm,lmaxFinal,lmaxTray,lim}; %#ok<AGROW>
        fprintf('h max estable (%s, %s): %.4g | limite teorico: %.4g\n', ...
            metH{im,1},metH{im,2},hm,lim);
    end
end

vn={'Epsilon','T','h','Metodo','Politica','Solver','Estado','ErrorMaxRef', ...
    'ErrorConservacion','Minimo','PasosAceptados','Rechazados','hEfectivo', ...
    'hMinAceptado','EvalF','Jacobianas','IterNewton','IterGS','RechGS','RechPos', ...
    'RechNewton','RechLinea','RhoMax','Segundos','RazonEspectral'};
tabla=cell2table(filas,'VariableNames',vn);
writetable(tabla,fullfile(carpetaRes,'MN_Grupo8_resultados_rigidez.csv'));
tablaH=cell2table(filasH,'VariableNames',{'Epsilon','Metodo','Solver', ...
    'hMaxEstable','LambdaMaxFinal','LambdaMaxTrayectoria','LimiteTeorico'});
writetable(tablaH,fullfile(carpetaRes,'MN_Grupo8_resultados_hmax_estable.csv'));
fprintf('\nGuardados: resultados_rigidez.csv y resultados_hmax_estable.csv\n');

%% 3. rho(T_GS) frente a h (condicion de convergencia de Gauss-Seidel)
q=p; q.eps=1e-2; refQ=ode15s(@(t,x)modelo(x,q),[0 100],x0,optRef);
estadosRho={x0, deval(refQ,50)}; etiqRho={'x_0 (inicial)','x(t=50)'};
hg=logspace(log10(0.1),log10(64),30);
filasR={};
for ig=1:numel(hg)
    fila=cell(1,5); fila{1}=hg(ig);
    for k=1:2
        W=eye(6)-hg(ig)*jacobiana_f(estadosRho{k},q);
        rho=max(abs(eig(-(tril(W)\triu(W,1)))));
        dd=all(abs(diag(W))>sum(abs(W-diag(diag(W))),2));   % dominancia diagonal estricta
        fila{2*k}=rho; fila{2*k+1}=double(dd);
    end
    filasR(end+1,:)=fila; %#ok<AGROW>
end
tablaR=cell2table(filasR,'VariableNames',{'h','Rho_x0','DomDiag_x0','Rho_x50','DomDiag_x50'});
writetable(tablaR,fullfile(carpetaRes,'MN_Grupo8_resultados_rho_gs.csv'));

%% 4. Figuras
% Figura 02: error frente al paso nominal; fila 1 sin reduccion, fila 2 con reduccion.
cfgNom={'Euler explicito','nominal','-','Euler explícito';
        'RK4','nominal','-','RK4';
        'Euler atras','nominal','GS','Euler atrás (GS)';
        'Euler atras','nominal','directo','Euler atrás (directo)'};
cfgRed={'Euler explicito','reduccion','-','Euler explícito';
        'RK4','reduccion','-','RK4';
        'RKF45','adaptativo','-','RKF45';
        'Euler atras','reduccion','GS','Euler atrás (GS)';
        'Euler atras','reduccion','directo','Euler atrás (directo)'};
colores=lines(5);
nomMet={'Euler explícito','RK4','RKF45','Euler atrás (GS)','Euler atrás (directo)'};  % color fijo por método
fig2=figure('Color','w','Name','Rigidez: error vs paso','Position',[40 40 1350 680]);
for ie=1:numel(epsilons)
    for fila=1:2
        ax=subplot(2,3,(fila-1)*3+ie); hold(ax,'on'); grid(ax,'on');
        if fila==1, C=cfgNom; tit='sin reducción de paso'; else, C=cfgRed; tit='con reducción de paso'; end
        hl=gobjects(0); et={};
        for c=1:size(C,1)
            m=strcmp(tabla.Metodo,C{c,1}) & strcmp(tabla.Politica,C{c,2}) & ...
              strcmp(tabla.Solver,C{c,3}) & tabla.Epsilon==epsilons(ie);
            sub=tabla(m,:); ok=strncmp(sub.Estado,'OK',2);
            ci=find(strcmp(nomMet,C{c,4}));
            hl(end+1)=plot(sub.h(ok),max(sub.ErrorMaxRef(ok),1e-12),'-o', ...
                'Color',colores(ci,:),'LineWidth',1.3,'MarkerSize',6); %#ok<AGROW>
            et{end+1}=C{c,4}; %#ok<AGROW>
            if any(~ok)
                plot(sub.h(~ok),0.6*0.45^(c-1)*ones(sum(~ok),1),'x','Color',colores(ci,:), ...
                    'MarkerSize',11,'LineWidth',1.8);   % altura distinta por método para que no se superpongan
            end
        end
        hx=plot(NaN,NaN,'kx','MarkerSize',10,'LineWidth',1.6);
        set(ax,'XScale','log','YScale','log'); if fila==1, ylim([1e-11 1]); else, ylim([1e-11 1e3]); end   % fila 2: margen superior para la leyenda
        xlim([0.4 5]);
        xticks(pasos); xticklabels(string(pasos)); set(ax,'XMinorTick','off');
        xlabel('Paso nominal h');
        if ie==1, ylabel('Error máximo frente a ode15s'); end
        title(sprintf('\\epsilon=%g, %s',epsilons(ie),tit));
        if ie==3
            if fila==1, legend([hl hx],[et {'inadmisible / fallo'}],'Location','southeast','FontSize',7);
            else, legend([hl hx],[et {'inadmisible / fallo'}],'Location','north','NumColumns',3,'FontSize',7); end
        end
    end
end
exportgraphics(fig2,fullfile(carpetaFiguras,'fig02_rigidez_estabilidad.png'),'Resolution',300);
fprintf('Figura exportada: fig02_rigidez_estabilidad.png\n');

% Figura 09: paso aceptado a lo largo de la integracion (epsilon=0.01, h nominal=4)
fig9=figure('Color','w','Name','Paso efectivo'); hold on; grid on;
col9=lines(numel(trayectorias)); hl9=gobjects(0); et9={};
for k=1:numel(trayectorias)
    tt=trayectorias{k}{3}; hh=diff(tt);
    if numel(tt)<2, continue, end
    hl9(end+1)=stairs(tt(1:end-1),hh,'LineWidth',1.4,'Color',col9(k,:)); %#ok<AGROW>
    nom=trayectorias{k}{1}; if ~strcmp(trayectorias{k}{2},'-'), nom=[nom ' (' trayectorias{k}{2} ')']; end
    et9{end+1}=strrep(nom,'atras','atrás'); %#ok<AGROW>
end
hn=yline(4,'--k','h nominal');
xlabel('Tiempo rápido t'); ylabel('Paso aceptado h');
title('Paso aceptado durante la integración (\epsilon=0.01, h nominal=4)');
legend([hl9 hn],[et9 {'h nominal = 4'}],'Location','best','FontSize',8);
% fig09 ya no se usa en el paper (se eliminó el Anexo C); no se exporta.
% exportgraphics(fig9,fullfile(carpetaFiguras,'fig09_paso_efectivo.png'),'Resolution',300);

% Figura 08: radio espectral de Gauss-Seidel frente a h
fig8=figure('Color','w','Name','rho GS vs h');
semilogx(tablaR.h,tablaR.Rho_x0,'-o','LineWidth',1.4); hold on;
semilogx(tablaR.h,tablaR.Rho_x50,'-s','LineWidth',1.4);
yline(1,'--k','\rho=1'); grid on;
set(gca,'YScale','log');            % escala log: se ve el cruce de rho=1 y el pico a la vez
xlabel('Paso h'); ylabel('\rho(T_{GS})');
title('Radio espectral de Gauss–Seidel para W = I - hJ_f (\epsilon=0.01)');
legend(etiqRho,'Location','northwest');
exportgraphics(fig8,fullfile(carpetaFiguras,'fig08_rho_gauss_seidel.png'),'Resolution',300);
fprintf('Figura exportada: fig08_rho_gauss_seidel.png\n');

fprintf('\nNota: con reduccion de paso comparar PasosAceptados y hEfectivo, no solo h nominal.\n');
fprintf('La razon espectral es un indicador local, no prueba por si sola rigidez.\n');

%% FUNCIONES LOCALES
function [fila,t]=ejecutar(metodo,politica,solver,x0,p,T,h,opt,ref,ratio)
    o=opt; o.reducir=strcmp(politica,'reduccion');
    if strcmp(solver,'directo'), o.solver='directo'; else, o.solver='gs'; end
    reloj=tic;
    switch metodo
        case 'Euler explicito', [t,X,info]=explicito(x0,p,T,h,1,o);
        case 'RK4',             [t,X,info]=explicito(x0,p,T,h,4,o);
        case 'RKF45',           [t,X,info]=rkf45(x0,p,T,h,h,o.rtol,o.atol,o);
        otherwise,              [t,X,info]=euler_atras(x0,p,T,h,o);
    end
    seg=toc(reloj); info=completar(info);
    if strcmp(info.estado,'OK')
        err=max(abs(X-deval(ref,min(max(t,0),T))),[],'all');
        cons=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
        mn=min(X,[],'all'); hEf=T/max(info.aceptados,1);
    else
        err=NaN; cons=NaN; mn=NaN; hEf=NaN;
    end
    if isinf(info.hMin), info.hMin=NaN; end
    fila={p.eps,T,h,metodo,politica,solver,info.estado,err,cons,mn,info.aceptados, ...
        info.rechazados,hEf,info.hMin,info.evalf,info.jac,info.newton,info.gs, ...
        info.rechGS,info.rechPos,info.rechNewton,info.rechLinea,info.rhoMax,seg,ratio};
end

function hm=biseccion_hmax(metodo,solver,x0,p,Tb,opt)
    o=opt; o.reducir=false; if strcmp(solver,'directo'), o.solver='directo'; else, o.solver='gs'; end
    lo=0.1; hi=64;
    if corre_nominal(metodo,x0,p,Tb,hi,o), hm=hi; return, end
    for k=1:14
        mid=sqrt(lo*hi);
        if corre_nominal(metodo,x0,p,Tb,mid,o), lo=mid; else, hi=mid; end
    end
    hm=lo;
end

function ok=corre_nominal(metodo,x0,p,T,h,o)
    switch metodo
        case 'Euler explicito', [~,~,info]=explicito(x0,p,T,h,1,o);
        case 'RK4',             [~,~,info]=explicito(x0,p,T,h,4,o);
        otherwise,              [~,~,info]=euler_atras(x0,p,T,h,o);
    end
    ok=strcmp(info.estado,'OK');
end

function info=completar(info)
    campos={'jac','newton','gs','rechGS','rechPos','rechNewton','rechLinea','rhoMax','evalf','hMin'};
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
