clear all
close all
clc


% Parámetros
kc = 3.398e-7;              % Coeficiente de transferencia [1/s]
a = 2.307e4;                % Superficie específica [1/m]
K = 0.9;

% Condiciones iniciales
Ca_f0 = 10.71;            % Concentración inicial en la fruta [g/L]
Ca_s0 = 0;            % Concentración inicial en el solvente [g/L]
Ca0 = [Ca_f0; Ca_s0];

% Tiempo de simulación
tspan = [0 3600];       % 1 hora

% Resolver con ode45
[t, Ca] = ode45(@(t, Ca) modelo_antocianina(t, Ca, kc, K, a), tspan, Ca0);

% Extraer concentraciones
Ca_f = Ca(:,1);
Ca_s = Ca(:,2);

% Calcular rendimiento

X = (Ca_f(1)-Ca_f)/Ca_f(1);

% Graficar resultados
figure(1)
hold on
plot(t/60, Ca_s, 'b', 'LineWidth', 2);

plot(t/60, Ca_f, 'r', 'LineWidth', 2);
xlabel('Tiempo [min]');
ylabel('Concentración [g/L]');
legend('Solvente','Fruta');
title('Transferencia de Antocianinas con Ultrasonido');
grid on;
hold off

figure(2)
hold on
plot(t/60, X,'b','LineWidth',2);
xlabel('Tiempo [min]');
ylabel('Rendimiento');
ylim([0 1])
xlim([0 10])
title('Caso mejorado Antocianinas: Rendimiento');
hold off
