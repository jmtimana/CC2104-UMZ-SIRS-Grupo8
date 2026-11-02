%% Grupo 8 - UMZ-SIRS | 09. Validacion de las soluciones de referencia (NUEVO)
% Responde a la observacion del docente: explicar como se valida la referencia en
% los casos de rigidez. Para 4 casos se comprueba:
%   (i)   convergencia de ode15s al estrechar RelTol (1e-8 ... 1e-12): diferencia
%         maxima contra la siguiente tolerancia y contra la mas estricta;
%   (ii)  contraste con integradores independientes (ode45 y ode113, no rigidos);
%   (iii) conservacion poblacional y positividad de la referencia;
%   (iv)  en los escenarios I y II, comparacion con el equilibrio ANALITICO.
% Salidas: ../resultados/MN_Grupo8_resultados_validacion_referencia.csv
%          ../resultados/MN_Grupo8_resultados_validacion_equilibrio.csv
clear; clc; close all; format long g

carpetaRes='../resultados';
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Casos
base.b1=1.1; base.b2=0.92; base.mu1=1; base.K=0.9;
base.beta=6; base.gamma=1.2; base.eta=0.7; base.mu2=1;
xAlt=[0.7;0.2;0.1;0.8;0.2;0];

casos(1).nombre='Base eps=1e-2 (T=100)'; casos(1).p=base; casos(1).p.eps=1e-2;
casos(1).T=100; casos(1).x0=xAlt; casos(1).escenario=0;
casos(2).nombre='Base eps=1e-3 (T=1000)'; casos(2).p=base; casos(2).p.eps=1e-3;
casos(2).T=1000; casos(2).x0=xAlt; casos(2).escenario=0;

bul.mu1=5.2143; bul.mu2=5.2143; bul.gamma=0.1667; bul.eta=0.1; bul.K=0.9; bul.eps=1/50;
I0=0.02;
casos(3).nombre='Escenario I (tau=30)'; casos(3).p=bul;
casos(3).p.b1=6.7041; casos(3).p.b2=5.9041; casos(3).p.beta=4.8429;
casos(3).T=1500; casos(3).x0=[0.8;0.2;0;1-I0;I0;0]; casos(3).escenario=1;
casos(4).nombre='Escenario II (tau=30)'; casos(4).p=bul;
casos(4).p.b1=5.9041; casos(4).p.b2=6.7041; casos(4).p.beta=5.919;
casos(4).T=1500; casos(4).x0=[0.8;0;0.2;1-I0;I0;0]; casos(4).escenario=2;

rtols=[1e-8 1e-9 1e-10 1e-11 1e-12];
nGrid=4001;
filas={}; filasEq={};

%% 2. Validacion por caso
for c=1:numel(casos)
    p=casos(c).p; T=casos(c).T; x0=casos(c).x0;
    tg=linspace(0,T,nGrid);
    fun=@(t,x) modelo(x,p);
    fprintf('\n===== %s =====\n',casos(c).nombre);

    % (i) ode15s con tolerancias decrecientes
    Y=cell(1,numel(rtols)); sols=cell(1,numel(rtols)); seg=zeros(1,numel(rtols));
    for k=1:numel(rtols)
        reloj=tic;
        sols{k}=ode15s(fun,[0 T],x0,odeset('RelTol',rtols(k),'AbsTol',1e-2*rtols(k)));
        seg(k)=toc(reloj); Y{k}=deval(sols{k},tg);
    end
    Yref=Y{end};
    for k=1:numel(rtols)
        errEstricta=max(abs(Y{k}-Yref),[],'all');
        if k<numel(rtols), errSig=max(abs(Y{k}-Y{k+1}),[],'all'); else, errSig=NaN; end
        cons=max([max(abs(sum(Y{k}(1:3,:),1)-1)),max(abs(sum(Y{k}(4:6,:),1)-1))]);
        [ns,nf]=estadisticas(sols{k});
        filas(end+1,:)={casos(c).nombre,'ode15s',rtols(k),1e-2*rtols(k),errEstricta,errSig, ...
            cons,min(Y{k},[],'all'),ns,nf,seg(k)}; %#ok<AGROW>
        fprintf('ode15s RelTol=%.0e | dif. vs estricta=%.3e | dif. vs siguiente=%.3e | cons=%.2e\n', ...
            rtols(k),errEstricta,errSig,cons);
    end

    % (ii) integradores independientes (error frente a ode15s mas estricta)
    ind={'ode45',@ode45; 'ode113',@ode113};
    for k=1:size(ind,1)
        try
            reloj=tic;
            sol=ind{k,2}(fun,[0 T],x0,odeset('RelTol',1e-12,'AbsTol',1e-14));
            s=toc(reloj); Yk=deval(sol,tg);
            [ns,nf]=estadisticas(sol);
            cons=max([max(abs(sum(Yk(1:3,:),1)-1)),max(abs(sum(Yk(4:6,:),1)-1))]);
            err=max(abs(Yk-Yref),[],'all');
            filas(end+1,:)={casos(c).nombre,ind{k,1},1e-12,1e-14,err,NaN,cons, ...
                min(Yk,[],'all'),ns,nf,s}; %#ok<AGROW>
            fprintf('%s RelTol=1e-12 | dif. vs ode15s estricta=%.3e\n',ind{k,1},err);
        catch ME
            filas(end+1,:)={casos(c).nombre,ind{k,1},1e-12,1e-14,NaN,NaN,NaN,NaN,NaN,NaN,NaN}; %#ok<AGROW>
            fprintf('%s fallo: %s\n',ind{k,1},ME.message);
        end
    end

    % (iv) equilibrio analitico (solo escenarios I y II, tau=30 ya cerca del equilibrio)
    if casos(c).escenario>0
        eq=equilibrio_analitico(p);
        if casos(c).escenario==1, e=eq.M; rama='M'; else, e=eq.Z; rama='Z'; end
        xf=Yref(:,end);
        vars={'U','M','Z','I'}; val=[xf(1) xf(2) xf(3) xf(5)];
        an=[e.U e.M e.Z e.I];
        for v=1:4
            filasEq(end+1,:)={casos(c).nombre,rama,vars{v},val(v),an(v),val(v)-an(v),e.ok}; %#ok<AGROW>
            fprintf('Equilibrio rama %s: %s numerico=%.8f analitico=%.8f dif=%.2e\n', ...
                rama,vars{v},val(v),an(v),val(v)-an(v));
        end
    end
end

%% 3. Exportar
tabla=cell2table(filas,'VariableNames',{'Caso','Integrador','RelTol','AbsTol', ...
    'Dif_vs_Estricta','Dif_vs_Siguiente','Error_Conservacion','Minimo','Pasos', ...
    'Evaluaciones','Segundos'});
writetable(tabla,fullfile(carpetaRes,'MN_Grupo8_resultados_validacion_referencia.csv'));
tablaEq=cell2table(filasEq,'VariableNames',{'Caso','Rama','Variable','Valor_ode15s', ...
    'Valor_Analitico','Diferencia','Rama_Coherente'});
writetable(tablaEq,fullfile(carpetaRes,'MN_Grupo8_resultados_validacion_equilibrio.csv'));
fprintf('\nGuardados: validacion_referencia.csv y validacion_equilibrio.csv\n');
fprintf('La referencia es numerica: Dif_vs_Siguiente cercana a Dif_vs_Estricta indica convergencia.\n');

%% FUNCIONES LOCALES
function [ns,nf]=estadisticas(sol)
    ns=NaN; nf=NaN;
    if isfield(sol,'stats')
        if isfield(sol.stats,'nsteps'), ns=sol.stats.nsteps; end
        if isfield(sol.stats,'nfevals'), nf=sol.stats.nfevals; end
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

function eq=equilibrio_analitico(p)
    Iend=@(be) max(0,(1-(p.gamma+p.mu2)./be)./(1+p.gamma/(p.eta+p.mu2)));
    Ms=max(0,1-p.mu1/p.b1); I_M=Iend(p.beta*(1+Ms));
    eq.M=struct('I',I_M,'U',1-Ms,'M',Ms,'Z',0, ...
        'ok',(p.b1>p.mu1)&&(p.b2/(1-p.K*I_M)<=p.b1));
    Zf=@(I) max(0,1-p.mu1*(1-p.K*I)/p.b2);
    g=@(I) I-Iend(p.beta/(1+Zf(I)));
    lo=0; hi=1;
    for k=1:60
        mid=(lo+hi)/2;
        if g(mid)<0, lo=mid; else, hi=mid; end
    end
    I_Z=(lo+hi)/2; Zz=Zf(I_Z);
    eq.Z=struct('I',I_Z,'U',1-Zz,'M',0,'Z',Zz, ...
        'ok',(p.b2/(1-p.K*I_Z)>p.mu1)&&(p.b1*(1-p.K*I_Z)/p.b2<=1));
end
