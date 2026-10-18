%% Grupo 8 - Modelo UMZ-SIRS: formulacion y validaciones basicas
% Metodos Numericos CC2104 - UTEC, 2026-2.
% Preparacion del avance 3. Ejecutar todas las secciones de arriba abajo.
% Este archivo NO integra la solucion en el tiempo; valida la funcion f(x).
% Para convertirlo a Live Script: abrir en MATLAB y guardar como .mlx.

clear; clc; close all;
format long g

%% 1. Parametros exploratorios (propuestos en el avance del grupo)
% Estos valores NO deben atribuirse al articulo original sin verificar.
p.b1 = 1.1;    % U -> M
p.b2 = 0.92;   % U -> Z
p.mu1 = 1.0;  % reemplazo/perdida de interes en la capa informativa
p.K = 0.9;    % sensibilidad de la tasa U -> Z a I (0 <= K < 1)
p.beta = 6.0; % contagio epidemiologico
p.gamma = 1.2;% recuperacion
p.eta = 0.7;  % perdida de inmunidad
p.mu2 = 1.0;  % demografia de SIRS
p.eps = 1e-2; % tau = eps*t, escala lenta

assert(p.K >= 0 && p.K < 1, 'K debe estar en [0,1).');
assert(p.eps > 0 && p.eps < 1, 'Elegir 0 < eps < 1.');

%% 2. Condiciones iniciales normalizadas
% Orden de las variables: x = [U; M; Z; S; I; R].
x0 = [0.8; 0.2; 0.0; 0.8; 0.2; 0.0];
x0_conZ = [0.7; 0.2; 0.1; 0.8; 0.2; 0.0];

assert(all(x0 >= 0) && all(x0_conZ >= 0), 'Poblaciones negativas.');
assert(abs(sum(x0(1:3))-1) < 1e-12, 'UMZ no suma 1.');
assert(abs(sum(x0(4:6))-1) < 1e-12, 'SIRS no suma 1.');
assert(abs(sum(x0_conZ(1:3))-1) < 1e-12, 'UMZ alternativo no suma 1.');
assert(abs(sum(x0_conZ(4:6))-1) < 1e-12, 'SIRS alternativo no suma 1.');

%% 3. Tasas de acoplamiento y derivadas iniciales
b2_efectiva = p.b2/(1-p.K*x0(5));
beta_efectiva = p.beta*(1+x0(2))/(1+x0(3));
f0 = modelo_umz_sirs(x0,p);

fprintf('b2 efectiva inicial: %.12g\n', b2_efectiva);
fprintf('beta efectiva inicial: %.12g\n', beta_efectiva);
variables = {'U'; 'M'; 'Z'; 'S'; 'I'; 'R'};
disp(table(variables, x0, f0, 'VariableNames', {'Variable','ValorInicial','DerivadaInicial'}));

%% 4. Verificacion de conservacion poblacional
% Cuando cada suma vale 1, la suma de derivadas debe valer 0.
cons_umz = sum(f0(1:3));
cons_sirs = sum(f0(4:6));
assert(abs(cons_umz) < 1e-12, 'Falla conservacion UMZ.');
assert(abs(cons_sirs) < 1e-12, 'Falla conservacion SIRS.');

% Prueba adicional: las sumas de derivadas deben satisfacer las identidades
% incluso para estados NO normalizados, mientras las tasas esten definidas.
x_prueba = [0.5; 0.2; 0.1; 0.7; 0.1; 0.1];
f_prueba = modelo_umz_sirs(x_prueba,p);
ident_umz = p.mu1*(1-sum(x_prueba(1:3)));
ident_sirs = p.eps*p.mu2*(1-sum(x_prueba(4:6)));
assert(abs(sum(f_prueba(1:3))-ident_umz) < 1e-12);
assert(abs(sum(f_prueba(4:6))-ident_sirs) < 1e-12);
fprintf('Conservacion: OK (errores iniciales %.3e y %.3e).\n', ...
    abs(cons_umz), abs(cons_sirs));

%% 5. Condicion necesaria de positividad del modelo continuo
% En cada frontera xi=0 se exige fi(x) >= 0. Se verifica para el
% conjunto de prueba; NO es una demostracion numerica universal.
tol_pos = 1e-10;
for i = 1:6
    x_frontera = x0_conZ;
    if i <= 3
        x_frontera(1:3) = [0; 0; 0];
        otros = setdiff(1:3,i);
        x_frontera(otros) = 0.5;
    else
        x_frontera(4:6) = [0; 0; 0];
        otros = setdiff(4:6,i);
        x_frontera(otros) = 0.5;
    end
    derivada = modelo_umz_sirs(x_frontera,p);
    assert(derivada(i) >= -tol_pos, 'Falla test de frontera en variable %d.',i);
end
fprintf('Pruebas de frontera de positividad: OK.\n');

%% 6. Advertencia sobre Z(0)=0
% dZ/dt = (b2/(1-K*I)*U - mu1)*Z, luego Z=0 es invariante.
% Para estudiar el efecto del escepticismo, usar x0_conZ.
f0_conZ = modelo_umz_sirs(x0_conZ,p);
disp(table(variables, x0_conZ, f0_conZ, ...
    'VariableNames', {'Variable','InicialConZ','DerivadaConZ'}));

%% 7. Resumen de validacion
fprintf('\nVALIDACION INICIAL COMPLETADA.\n');
fprintf('No se ha realizado integracion temporal ni comparacion de metodos.\n');
fprintf('Siguiente etapa: Euler explicito, RK4, RKF45 y Euler hacia atras.\n');

%% Funcion local: lado derecho del sistema en el tiempo rapido t
function dx = modelo_umz_sirs(x,p)
    assert(numel(x)==6, 'x debe contener 6 componentes.');
    x=x(:);
    U=x(1); M=x(2); Z=x(3);
    S=x(4); I=x(5); R=x(6);
    assert(1-p.K*I > 0, 'Denominador 1-K*I no positivo.');
    assert(1+Z > 0, 'Denominador 1+Z no positivo.');

    b2_bar = p.b2/(1-p.K*I);
    beta_bar = p.beta*(1+M)/(1+Z);
    dx = zeros(6,1);
    dx(1) = p.mu1 - p.b1*U*M - b2_bar*U*Z - p.mu1*U;
    dx(2) = p.b1*U*M - p.mu1*M;
    dx(3) = b2_bar*U*Z - p.mu1*Z;
    dx(4) = p.eps*(p.mu2-beta_bar*S*I+p.eta*R-p.mu2*S);
    dx(5) = p.eps*(beta_bar*S*I-(p.gamma+p.mu2)*I);
    dx(6) = p.eps*(p.gamma*I-(p.eta+p.mu2)*R);
end
