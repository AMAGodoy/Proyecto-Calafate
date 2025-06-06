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

% Parámetros
Kc        = 6.0793e-7;  % [m/s]      Coeficiente de transferencia de masa global
a         = 1.5e3;      % [m^2/m^3]  Razón de superficie de hollejo por volumen de reactor
K         = 0.2;        % [-]        Coeficiente de partición de antocianinas en el hollejo
v         = 1.15e-4;    % [m/s]      Velocidad promedio del flujo
r_reactor = 0.0175;     % [m]        Radio basal del reactor
L_reactor = 0.1;        % [m]        Largo del reactor
Ca0       = 10.71;      % [kg/m^3]   Concentracion inicial de antocianina en un hollejo

%% Vector Eje Z
z_0 = 0;                % [m]
dz  = L_reactor/60;     % [m]
z_f = L_reactor;        % [m]
Z   = [z_0:dz:z_f];

%% Vector Eje R
r_0 = 0;                % [m]
dr = r_reactor/20;      %[m]
r_f = 2*r_reactor;      %[m]
R = [r_0+dr:dr:r_f+dr]; % Uso de un valor en r al incio y al final para ver paredes del tubo

%% Vector Temporal
t_0 = 0;                % [s]
dt  = 60;               % [s]
t_f = 6*3600;           % [s]
T   = [t_0:dt:t_f];



% Inicio de matrices
%% Ca Solvente
Ca_S = zeros(length(Z), length(R), length(T));

%% Ca de un fruto
Ca_F = zeros(length(Z), length(R), length(T));

%% Ca superficie fruto
Ca_SF = zeros(length(Z), length(R), length(T));

%% Flux entre superficie de fruto y solvente
Ja = zeros(length(Z), length(R), length(T));

%% Rendimiento
X  = zeros(length(T),1);



% Condiciones iniciales
%% Concentracion en el solvente
Ca_S(1:end,2:end-1,1) = 0; % [kg/m^3]

%% Concentracion en el fruto
Ca_F(1:end,2:end-1,1) = Ca0; % [kg/m^3]

%% Concentracion en la superficie del fruto
Ca_SF(1:end,2:end-1,1) = Ca_F(1:end,2:end-1,1) * K; % [kg/m^3]

%% Flux desde la superficie al solvente
Ja(1:end,1) = 0; % [kg/(s*m^2)]



% Diferencias finitas
%% Constantes
c1 = 1/(1 + v*dt/dz + Kc*a*dt); % [-]
c2 = v/(dz/dt + v + Kc*a*dz) ;   % [-]
c3 = Kc*a/(1/dt+v/dz+Kc*a);     % [-]

%% Evaluacion de las diferencias
for k = 2:length(T)
    %% Cálculo de valores de concentracion en solvente en nuevo tiempo
    for i = 2:length(Z)-1
        Ja(i,2:end-1,k) = Kc * (Ca_SF(i,2:end-1,k-1) - Ca_S(i,2:end-1,k-1));
        Ca_F(i,2:end-1,k) = Ca_F(i,2:end-1,k-1) - Ja(i,2:end-1,k) * a * dt;
        Ca_SF(i,2:end-1,k) = K * Ca_F(i,2:end-1,k);

        Ca_S(i,2:end-1,k) = Ca_S(i,2:end-1,k-1)*c1 + Ca_S(i-1,2:end-1,k)*c2 + Ca_SF(i,2:end-1,k)*c3;
    end
end



% Calculo de rendimiento
M_tot = sum(sum(Ca_F));

X = (M_tot(1)-M_tot)./M_tot(1);



% Graficar
ZZ = -Z;

%% Graficar concentracion en el solvente
figure(1)
hold on
for k = 1:length(T)
    t_hora = int2str(floor(T(k)/3600));
    t_minuto = int2str(mod(T(k)/60,60));
    imagesc(R,ZZ,Ca_S(:,:,k))
    xlabel("Radio [m]")
    ylabel("Largo [m]")
    h = colorbar();
    xlim("tight")
    ylim("tight")
    title(cstrcat("Caso base Antocianinas (Medio Solvente), ",t_hora, ' [h] y ', t_minuto, ' [min]'))
    title(h,"Concentración [kg/m^3]")
    pause(0.01)
end
hold off

%% Graficar concentracion en el fruto
figure(2)
hold on
for k = 1:length(T)
    t_hora = int2str(floor(T(k)/3600));
    t_minuto = int2str(mod(T(k)/60,60));
    imagesc(R,ZZ,Ca_F(:,:,k))
    xlabel("Radio [m]")
    ylabel("Largo [m]")
    h = colorbar();
    xlim("tight")
    ylim("tight")
    title(cstrcat("Caso base Antocianinas (Medio Fruto): ",t_hora, ' [h] y ', t_minuto, ' [min]'))
    title(h,"Concentración [kg/m^3]")
    pause(0.01)
end
hold off

%% Graficar rencimiento en el tiempo
figure(3)
hold on
plot(T/60,X,'b','LineWidth',2)
xlabel("Tiempo [min]")
ylabel("Rendimiento")
xlim("tight")
title("Caso base Antocianinas: Rendimiento")
grid on
hold off
