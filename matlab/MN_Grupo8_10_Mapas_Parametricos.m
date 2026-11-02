%% Grupo 8 - UMZ-SIRS | 10. Mapas parametricos 2D (NUEVO)
% Objetivo especifico 4: condiciones parametricas bajo las cuales la difusion de
% desinformacion transforma una epidemia en extincion en una endemia persistente.
% Mapa A: b1 x beta (b2 base).   Mapa B: b1 x b2 con beta fijo cuyo R00<1.
% En cada punto: I final por RK4 (vectorizado sobre todos los puntos), I* ANALITICO,
% regimen informativo y clasificacion extincion/endemia (numerica y analitica),
% contrastada con el control sin informacion (I* de beta).
% Salidas: ../resultados/MN_Grupo8_resultados_mapas.csv
%          ../figuras/fig06_mapas_parametricos.png
% Nota: guardar en UTF-8 (los textos de las figuras llevan tildes).
clear; clc; close all; format long g

carpetaFiguras='../figuras'; carpetaRes='../resultados';
if ~exist(carpetaFiguras,'dir'), mkdir(carpetaFiguras); end
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

%% 1. Parametros y rejilla
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=0.01;
x0=[0.7;0.2;0.1;0.8;0.2;0];
N=31;                          % puntos por eje (31x31 por mapa)
tauFinal=20; T=tauFinal/p.eps; h=0.25;
umbral=1e-3;                   % I<umbral se clasifica como extincion
betaMapaB=2;                   % R00=beta/(gamma+mu2)=0.909<1 en el conjunto base

b1v=linspace(0.5,3,N);
betav=linspace(1,8,N);
b2v=linspace(0.25,3,N);
mapas(1).nombre='A: b1 x beta'; mapas(1).x=b1v; mapas(1).y=betav;
mapas(1).xl='b_1'; mapas(1).yl='\beta';
mapas(2).nombre='B: b1 x b2'; mapas(2).x=b1v; mapas(2).y=b2v;
mapas(2).xl='b_1'; mapas(2).yl='b_2';

Iend=@(be) max(0,(1-(p.gamma+p.mu2)./be)./(1+p.gamma/(p.eta+p.mu2)));
filas={}; resMapa=cell(1,2);

%% 2. Calculo
for m=1:2
    [X1,Y1]=meshgrid(mapas(m).x,mapas(m).y);   % filas: eje y; columnas: eje x
    npts=numel(X1);
    B1=X1(:).'; 
    if m==1, B2=p.b2*ones(1,npts); BE=Y1(:).'; else, B2=Y1(:).'; BE=betaMapaB*ones(1,npts); end
    [Ifin,estado]=rk4_vectorizado(x0,p,B1,B2,BE,T,h);
    Ian=NaN(1,npts); reg=cell(1,npts); Ictl=NaN(1,npts);
    for k=1:npts
        q=p; q.b1=B1(k); q.b2=B2(k); q.beta=BE(k);
        eq=equilibrio_analitico(q);
        Ian(k)=eq.I; reg{k}=eq.regimen; Ictl(k)=Iend(BE(k));
    end
    claseNum=Ifin>=umbral; claseAn=Ian>=umbral; claseCtl=Ictl>=umbral;
    for k=1:npts
        if ~estado(k), cn='INADMISIBLE'; elseif claseNum(k), cn='endemia'; else, cn='extincion'; end
        if isnan(Ian(k)), ca='indefinido'; elseif claseAn(k), ca='endemia'; else, ca='extincion'; end
        if claseCtl(k), cc='endemia'; else, cc='extincion'; end
        if ~claseCtl(k) && strcmp(ca,'endemia'), efecto='endemia inducida por informacion';
        elseif claseCtl(k) && strcmp(ca,'extincion'), efecto='extincion inducida por informacion';
        else, efecto='sin cambio'; end
        concuerda=strcmp(cn,ca);
        filas(end+1,:)={mapas(m).nombre,B1(k),B2(k),BE(k),Ifin(k),Ian(k),Ictl(k), ...
            reg{k},cn,ca,cc,efecto,double(concuerda),double(estado(k))}; %#ok<AGROW>
    end
    resMapa{m}.I=reshape(Ifin,size(X1)); resMapa{m}.Ian=reshape(Ian,size(X1));
    resMapa{m}.X=X1; resMapa{m}.Y=Y1;
    resMapa{m}.ind=reshape(strcmp(reg,'indefinido'),size(X1));   % puntos sin rama pura coherente
    ok=estado & ~isnan(Ian);
    fprintf('Mapa %s: concordancia numerico/analitico = %.1f %% (puntos admisibles: %d/%d)\n', ...
        mapas(m).nombre,100*mean(claseNum(ok)==claseAn(ok)),nnz(estado),npts);
end

tabla=cell2table(filas,'VariableNames',{'Mapa','b1','b2','beta','I_Final','I_Equilibrio_Analitico', ...
    'I_Control_Analitico','Regimen_Informativo','Clase_Numerica','Clase_Analitica', ...
    'Clase_Control','Efecto_Informacion','Concuerda','Admisible'});
writetable(tabla,fullfile(carpetaRes,'MN_Grupo8_resultados_mapas.csv'));
fprintf('Guardado: MN_Grupo8_resultados_mapas.csv\n');

%% 3. Figura
fig=figure('Color','w','Name','Mapas parametricos','Position',[40 40 1250 500]);
tit={sprintf('(a) I final en \\tau=%g (b_2=%g)',tauFinal,p.b2), ...
     sprintf('(b) I final en \\tau=%g (\\beta=%g, R_{00}<1)',tauFinal,betaMapaB)};
for m=1:2
    subplot(1,2,m);
    imagesc(mapas(m).x,mapas(m).y,resMapa{m}.I); set(gca,'YDir','normal'); hold on;
    colorbar; colormap(parula);
    [~,hc]=contour(resMapa{m}.X,resMapa{m}.Y,resMapa{m}.Ian,[umbral umbral],'w','LineWidth',2);
    hi=plot(resMapa{m}.X(resMapa{m}.ind),resMapa{m}.Y(resMapa{m}.ind),'kx','MarkerSize',4);
    if m==1
        yline(p.gamma+p.mu2,'--k','R_{00}=1','LabelHorizontalAlignment','left','LineWidth',1.2);
    end
    xlabel(mapas(m).xl); ylabel(mapas(m).yl); title(tit{m});
    lg=legend([hc hi],{'umbral analítico I^*=10^{-3}','régimen indefinido (sin rama pura)'}, ...
        'Location','northwest','FontSize',8);
    lg.Color=[0.55 0.55 0.55]; lg.TextColor='w';   % fondo gris para que se vea la línea blanca
end
annotation('textbox',[0.30 0.005 0.4 0.05],'String', ...
    'Línea blanca: umbral analítico I^*=10^{-3}; cruces: puntos sin régimen analítico definido', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',9);
exportgraphics(fig,fullfile(carpetaFiguras,'fig06_mapas_parametricos.png'),'Resolution',300);
fprintf('Figura exportada: fig06_mapas_parametricos.png\n');
fprintf('NOTA: clasificacion de horizonte finito (tau=%g); cerca del umbral puede discrepar.\n',tauFinal);

%% FUNCIONES LOCALES
function [Ifin,ok]=rk4_vectorizado(x0,p,B1,B2,BE,T,h)
    npts=numel(B1); n=round(T/h);
    X=repmat(x0,1,npts); mn=min(X,[],1); maxKI=p.K*X(5,:); bad=false(1,npts);
    for k=1:n
        a=fv(X,p,B1,B2,BE); b=fv(X+h*a/2,p,B1,B2,BE);
        c=fv(X+h*b/2,p,B1,B2,BE); d=fv(X+h*c,p,B1,B2,BE);
        X=X+h*(a+2*b+2*c+d)/6;
        bad=bad | any(~isfinite(X),1);
        mn=min(mn,min(X,[],1)); maxKI=max(maxKI,p.K*X(5,:));
    end
    Ifin=X(5,:);
    ok=~bad & mn>=-1e-10 & maxKI<1;
end

function dX=fv(X,p,B1,B2,BE)
    U=X(1,:); M=X(2,:); Z=X(3,:); S=X(4,:); I=X(5,:); R=X(6,:);
    a=B2./(1-p.K*I); b=BE.*(1+M)./(1+Z);
    dX=[p.mu1-B1.*U.*M-a.*U.*Z-p.mu1*U;
        B1.*U.*M-p.mu1*M;
        a.*U.*Z-p.mu1*Z;
        p.eps*(p.mu2-b.*S.*I+p.eta*R-p.mu2*S);
        p.eps*(b.*S.*I-(p.gamma+p.mu2)*I);
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
    eq.I0=Iend(p.beta);
    if eq.M.ok, eq.regimen='M'; eq.I=eq.M.I;
    elseif eq.Z.ok, eq.regimen='Z'; eq.I=eq.Z.I;
    elseif p.b1<=p.mu1 && p.b2/(1-p.K*eq.I0)<=p.mu1, eq.regimen='ninguno'; eq.I=eq.I0;
    else, eq.regimen='indefinido'; eq.I=NaN;
    end
end
