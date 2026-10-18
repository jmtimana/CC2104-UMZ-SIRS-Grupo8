%% Grupo 8: analisis paramétrico exploratorio UMZ-SIRS
% PENDIENTE DE EJECUTAR EN MATLAB. No inventar resultados.
% Un parametro por vez; los demas se conservan en los valores base.
clear; clc; close all;
p.b1=1.1; p.b2=0.92; p.mu1=1; p.K=0.9;
p.beta=6; p.gamma=1.2; p.eta=0.7; p.mu2=1; p.eps=0.01;
x0=[0.7;0.2;0.1;0.8;0.2;0];
T=1000; h=0.1;  % tau final = 10; paso en tiempo rapido
n=round(T/h);
parametros={'b1','b2','K','beta'};
valores={[0.5 1.1 2.0],[0.5 0.92 2.0],[0 0.3 0.6 0.9],[2 4 6 8]};
filas={};
for j=1:numel(parametros)
    nombre=parametros{j};
    for v=valores{j}
        q=p; q.(nombre)=v;
        X=zeros(6,n+1); X(:,1)=x0;
        valido=true;
        for k=1:n
            x=X(:,k);
            if ~admisible(x,q), valido=false; break; end
            a=f(x,q); b=f(x+h*a/2,q); c=f(x+h*b/2,q); d=f(x+h*c,q);
            X(:,k+1)=x+h*(a+2*b+2*c+d)/6;
            if ~admisible(X(:,k+1),q), valido=false; break; end
        end
        if valido
            conserv=max([max(abs(sum(X(1:3,:),1)-1)),max(abs(sum(X(4:6,:),1)-1))]);
            ifinal=X(5,end); ipico=max(X(5,:));
            % 'Rebrote' no se determina solo con el pico; requiere inspeccion de I(t).
            fprintf('%s=%g: I_final=%g; I_pico=%g; conserv=%g\n',nombre,v,ifinal,ipico,conserv);
            filas(end+1,:)={nombre,v,ifinal,ipico,conserv,'OK'}; %#ok<AGROW>
        else
            filas(end+1,:)={nombre,v,NaN,NaN,NaN,'INADMISIBLE'}; %#ok<AGROW>
        end
    end
end
resultado=cell2table(filas,'VariableNames',{'Parametro','Valor','I_Final','I_Pico','Error_Conservacion','Estado'});
writetable(resultado,'MN_Grupo8_resultados_parametricos_PENDIENTE.csv');
disp(resultado);

function ok=admisible(x,p)
    ok=all(isfinite(x)) && all(x>=-1e-10) && 1-p.K*x(5)>0 && 1+x(3)>0;
end
function dx=f(x,p)
    U=x(1); M=x(2); Z=x(3); S=x(4); I=x(5); R=x(6);
    a=p.b2/(1-p.K*I); b=p.beta*(1+M)/(1+Z);
    dx=[p.mu1-p.b1*U*M-a*U*Z-p.mu1*U;
        p.b1*U*M-p.mu1*M;
        a*U*Z-p.mu1*Z;
        p.eps*(p.mu2-b*S*I+p.eta*R-p.mu2*S);
        p.eps*(b*S*I-(p.gamma+p.mu2)*I);
        p.eps*(p.gamma*I-(p.eta+p.mu2)*R)];
end
