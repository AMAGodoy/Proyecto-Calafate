clear all
close all
clc

% Parámetros
Kc = 3.398e-4;
a  = 2.307e4;
K  = 0.151;
v  = 1.083;
r  = 1.3e-4;
V_pelota = 4/3*pi()*r^3;
A_pelota = 4*pi()*r^2;
Ca0 = 10.71;

% Vector Z
z_0 = 0;
dz  = 0.5;
z_f = 20;
Z   = [z_0:dz:z_f];

% Vector R
r_0 = 0;
dr = r/10;
r_f = 2*r;
R = [r_0+dr:dr:r_f+dr]; % Uso de un valor en r al incio y al final para ver paredes del tubo -> cambio de coordenadas

% Vector T
t_0 = 0;
dt  = 0.5;
t_f = 40;
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
Ca_S(1:end,2:end-1,1) = 0;

%% Concentracion en el fruto
Ca_F(1:end,2:end-1,1) = Ca0;

%% Concentracion en la superficie del fruto
Ca_SF(1:end,2:end-1,1) = Ca_F(1:end,2:end-1,1) * K;

%% Flux desde la superficie al solvente
Ja(1:end,1) = 0;


% Diferencia
%% Constantes
c1 = 1/(1 + v*dt/dz + Kc*a*dt);
c2 = 1/(dz/dt + v + Kc*a*dz);
c3 = Kc*a/(1/dt+v/dz+Kc*a);

%c4 = 1/(1 + v*dt/dz);
%c5 = 1/(dz/dt + v);


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

figure(1)
hold on
for k = 1:length(T)
    t_tiempo = int2str(T(k));
    imagesc(R,ZZ,Ca_S(:,:,k))
    xlabel("Largo")
    ylabel("Radio")
    h = colorbar();
    xlim("tight")
    ylim("tight")
    title(cstrcat("Caso base Antocianinas, ",t_tiempo, " [s]"))
    title(h,"Concentración")
    pause(0.25)
end
hold off
