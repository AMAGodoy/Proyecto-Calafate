% Integrantes:
%			-Andrés Aguilera J.
%			-Nicolás Ebner Z.
%     -Ciro Fernandez F.
%     -Matías Godoy A.
%     -Nicolás Prokes H.
%     -Benjamín Ramos V.


clear all
close all
clc

% Parámetros (ver tablas y nomenclatura del paper)
ka_s = 1.36e-6;   % [m/s] coef. de transferencia fase solvente a 40C
ka_f = 1.36e-6;   % [m/s] coef. de transferencia fase fruta a 40C
a    = 1.5e3;     % [1/m] superficie específica de masa
e  = 0.95;        % [-] fracción volumétrica del solvente
K    = 0.325;     % [-] constante de partición a 45C

% Condiciones iniciales
Ca_s0 = 0.0;      % [kg/m3] concentración inicial en solvente
Ca_f0 = 10.71;    % [kg/m3] concentración inicial en fruta
Ca0   = [Ca_s0; Ca_f0];

% Tiempo de simulación: de 0 a 30 min (0-1800 s)
tspan = [0, 1800];

% Resolver con solver ode45
[t, Ca] = ode45(@(t, Ca) modelo_UAE(t, Ca, e, ka_s, ka_f, a, K), tspan, Ca0);

% Extraer resultados
Ca_s = Ca(:,1);
Ca_f = Ca(:,2);

% Calculo de rendimiento
rendimiento = ((Ca_f0 - Ca_f) / Ca_f0);

% Gráfica ode
figure(1)
hold on
plot(t/60, Ca_s, 'b-', 'LineWidth', 2)
plot(t/60, Ca_f, 'r-', 'LineWidth', 2)
legend('Ca_s (solvente)', 'Ca_f (fruta)')
xlabel('Tiempo [min]')
ylabel('Concentración [kg/m^3]')
title('Extracción de antocianinas con modelo de interfaz')
grid on
hold off

% Gráfica rendimiento
figure(2)
hold on
plot(t/60, rendimiento, 'g-', 'LineWidth', 2)
xlabel('Tiempo [min]')
ylabel('Rendimiento de extracción')
title('Rendimiento de extracción de antocianinas en fruta')
grid on
hold off

