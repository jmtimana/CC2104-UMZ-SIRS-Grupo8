%% Grupo 8: analisis parametrico UMZ-SIRS (ENTREGA FINAL)
% Un parametro por vez; los demas se conservan en los valores base.
% Integrador: RK4 de paso fijo (h=0.1). Horizonte: tau=20 (T=2000); se guarda
% tambien I en tau=10 para reproducir los valores del avance 3.
% Para cada caso se registra: I final, I en tau=10, pico y su tiempo, numero de
% rebrotes (aumentos de I >= delta tras haber disminuido), el equilibrio
% endemico ANALITICO (independiente de la integracion), el regimen informativo
% predicho (M, Z, ninguno) y el estado numerico.
% Salida: ../resultados/MN_Grupo8_resultados_parametricos.csv
clear; clc; close all; format long g

carpetaRes='../resultados';
if ~exist(carpetaRes,'dir'), mkdir(carpetaRes); end

p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=0.01;
x0=[0.7;0.2;0.1;0.8;0.2;0];
tauFinal=20; T=tauFinal/p.eps; h=0.1;
parametros={'b1','b2','K','beta'};
valores={[0.5 0.8 1.0 1.1 1.5 2.0 3.0], ...
         [0.25 0.5 0.92 1.5 2.0 3.0], ...
         [0 0.15 0.3 0.45 0.6 0.75 0.9], ...
         [1 2 3 4 5 6 7 8]};

filas={};
for j=1:numel(parametros)
    nombre=parametros{j};
    for v=valores{j}
        q=p; q.(nombre)=v;
        R00=q.beta/(q.gamma+q.mu2);
        eq=equilibrio_analitico(q);
        res=simular_rk4(x0,q,T,h);
        if res.valido
            I=res.I;
            ifinal=I(end); ipico=max(I);
            [~,kp]=max(I); tpico=(kp-1)*h*q.eps;
            n10=round(10/(q.eps*h))+1; itau10=I(min(n10,numel(I)));
            delta=max(1e-8,1e-3*ipico);
            nReb=contar_rebrotes(I,delta);
            difAn=ifinal-eq.I;
            fprintf('%s=%g: I_final=%.6g (analitico %.6g, regimen %s); I_pico=%.4g; rebrotes=%d\n', ...
                nombre,v,ifinal,eq.I,eq.regimen,ipico,nReb);
            filas(end+1,:)={nombre,v,ifinal,itau10,ipico,tpico,nReb,eq.I,difAn, ...
                eq.regimen,R00,res.cons,res.minimo,'OK'}; %#ok<AGROW>
        else
            filas(end+1,:)={nombre,v,NaN,NaN,NaN,NaN,NaN,eq.I,NaN,eq.regimen,R00, ...
                NaN,NaN,'INADMISIBLE'}; %#ok<AGROW>
        end
    end
end
resultado=cell2table(filas,'VariableNames',{'Parametro','Valor','I_Final','I_Tau10', ...
    'I_Pico','Tiempo_Pico_tau','Num_Rebrotes','I_Equilibrio_Analitico', ...
    'Dif_Final_vs_Analitico','Regimen_Informativo','R00_sin_informacion', ...
    'Error_Conservacion','Minimo','Estado'});
writetable(resultado,fullfile(carpetaRes,'MN_Grupo8_resultados_parametricos.csv'));
disp(resultado);
fprintf('\nGuardado: MN_Grupo8_resultados_parametricos.csv\n');
fprintf('Rebrote = aumento de I >= delta despues de una disminucion >= delta.\n');
fprintf('Si Dif_Final_vs_Analitico no es pequena, el horizonte tau=%g aun no alcanza el equilibrio.\n',tauFinal);

%% FUNCIONES LOCALES
function res=simular_rk4(x0,q,T,h)
    n=round(T/h); I=zeros(1,n+1); I(1)=x0(5); x=x0;
    res.valido=true; cons=0; mn=min(x0); nOk=1;
    for k=1:n
        a=modelo(x,q); b=modelo(x+h*a/2,q); c=modelo(x+h*b/2,q); d=modelo(x+h*c,q);
        x=x+h*(a+2*b+2*c+d)/6;
        if ~(all(isfinite(x)) && all(x>=-1e-10) && 1-q.K*x(5)>0 && 1+x(3)>0)
            res.valido=false; break
        end
        nOk=nOk+1; I(nOk)=x(5);
        cons=max(cons,max(abs(sum(x(1:3))-1),abs(sum(x(4:6))-1)));
        mn=min(mn,min(x));
    end
    res.I=I(1:nOk); res.cons=cons; res.minimo=mn;
end

function nReb=contar_rebrotes(I,delta)
    % Histeresis (zigzag): un rebrote es un aumento >= delta posterior a una
    % disminucion >= delta. Un ascenso monotono hacia el equilibrio no cuenta.
    dir=0; ref=I(1); nReb=0; hubo=false;
    for k=2:numel(I)
        if dir>=0
            if I(k)>ref, ref=I(k);
            elseif ref-I(k)>=delta, dir=-1; hubo=true; ref=I(k);
            end
        else
            if I(k)<ref, ref=I(k);
            elseif I(k)-ref>=delta
                dir=1; if hubo, nReb=nReb+1; end
                ref=I(k);
            end
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

function eq=equilibrio_analitico(p)
    % I*=(1-(gamma+mu2)/beta_ef)/(1+gamma/(eta+mu2)), I*>=0 (ver archivo 07).
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
