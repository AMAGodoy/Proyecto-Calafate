clear all
close all
clc

% Parámetros
Kc = 3.398e-8;  %[m/s]
a  = 2.307e4; %[m^2/m^3]
K  = 0.151; % [-]
v  = 0.0083; % [m/s]
r_reactor = 0.0175; % [m]
r_pelota  = 1.3e-4; % [m]
V_pelota = 4/3*pi()*r_pelota^3; % [m^3]
A_pelota = 4*pi()*r_pelota^2; % [m^2]
Ca0 = 10.71; % [kg/m^3]

% Vector Altura reactor
z_0 = 0; % [m]
dz  = 0.005; % [m]
z_f = 0.05; % [m]
Z   = [z_0:dz:z_f];

% Vector Radio reactor
r_0 = 0; % [m]
dr = r_reactor/20; %[m]
r_f = 2*r_reactor; %[m]
R = [r_0+dr:dr:r_f+dr]; % Uso de un valor en r al incio y al final para ver paredes del tubo -> cambio de coordenadas

% Vector T
t_0 = 0; % [s]
dt  = 60; % [s]
t_f = 7*3600; % [s]
T   = [t_0:dt:t_f];

% Inicio de matrices
%% Ca Solvente
Ca_S = zeros(length(Z), length(R), length(T));

%% Ca de un fruto
Ca_F = zeros(length(Z), length(R), length(T));

%% Ca superficie fruto
Ca_SF = zeros(length(Z), length(R), length(T));

% Flux entre superficie de fruto y solvente
Ja = zeros(length(Z), length(R), length(T));


Ca_SFtest = zeros(length(Z), length(R), length(T));



% Condición inicial
%% Concentracion en el solvente
Ca_S(1:end,2:end-1,1) = 0; % [kg/m^3]

%% Concentracion en el fruto
Ca_F(1:end,2:end-1,1) = Ca0; % [kg/m^3]

%% Concentracion en la superficie del fruto
Ca_SF(1:end,2:end-1,1) = Ca_F(1:end,2:end-1,1) * K; % [kg/m^3]

%% Flux desde la superficie al solvente
Ja(1:end,1) = 0; % [kg/(s*m^2)]


% Diferencia
%% Constantes
c1 = 1/(1 + v*dt/dz + Kc*a*dt); %  [-]
c2 = v/(dz/dt + v + Kc*a*dz); % [-]
c3 = Kc*a/(1/dt+v/dz+Kc*a); % [-]


% Diferencias finitas
for k = 2:length(T)
    %% Cálculo de valores de concentracion en solvente en nuevo tiempo
    for i = 2:length(Z)-1
        Ja(i,2:end-1,k) = Kc * (Ca_SF(i,2:end-1,k-1) - Ca_S(i,2:end-1,k-1));
        Ca_F(i,2:end-1,k) = Ca_F(i,2:end-1,k-1) - Ja(i,2:end-1,k) * A_pelota/V_pelota * dt;
        Ca_SF(i,2:end-1,k) = K * Ca_F(i,2:end-1,k);

        Ca_S(i,2:end-1,k) = Ca_S(i,2:end-1,k-1)*c1 + Ca_S(i-1,2:end-1,k)*c2 + Ca_SF(i,2:end-1,k)*c3;
    end
    %Ca_S(end,2:end-1,k) = Ca_S(end,2:end-1,k-1)*c4 + Ca_S(end-1,2:end-1,k)*c5;

    %% Cálculo de nuevo valor de concentraciones en el fruto en nuevo tiempo

end

% Graficar
ZZ = -Z;

%figure(1)
%hold on
%for k = 1:length(T)
%    t_tiempo = int2str(T(k));
%    imagesc(R,ZZ,Ca_S(:,:,k))
%    xlabel("Largo")
%    ylabel("Radio")
%    h = colorbar();
%    xlim("tight")
%    ylim("tight")
%    title(cstrcat("Caso base Antocianinas, ",t_tiempo, " [s]"))
%    title(h,"Concentración")
%    pause(0.25)
%end
%hold off


%figure(2)
%hold on
%for k = 1:length(T)
%    t_tiempo = int2str(T(k));
%    subplot(1,2,1)
%    imagesc(ZZ,R,Ca_S(:,:,k))
%    xlabel("Radio")
%    ylabel("Largo")
%    h = colorbar();
%    xlim("tight")
%    ylim("tight")
%   title(cstrcat("Caso base Antocianinas Solvente, ",t_tiempo, " [s]"))
%    title(h,"Concentración")

%    subplot(1,2,2)
%    imagesc(ZZ,R,Ca_F(:,:,k))
%    xlabel("Radio")
%    ylabel("Largo")
%   hi = colorbar();
%    xlim("tight")
%    ylim("tight")
%    title(cstrcat("Caso base Antocianinas Fruto, ",t_tiempo, " [s]"))
%    title(h,"Concentración")
%    pause(0.05)
%end
%hold off


figure(3)
hold on
for k = 1:length(T)
    t_hora = int2str(floor(T(k)/3600));
    t_minuto = int2str(mod(T(k)/60,60));
    t_tiempo = int2str(T(k));
    imagesc(R,ZZ,Ca_F(:,:,k))
    xlabel("Radio")
    ylabel("Largo")
    h = colorbar();
    xlim("tight")
    ylim("tight")
    title(cstrcat("Caso base Antocianinas (Medio Fruto): ",t_hora, ' [h] y ', t_minuto, ' [min]'))
    title(h,"Concentración")
    pause(0.05)
end
hold off
