close all
clear all
clc
%% Dados
H = 0.24; 
W = 0.24; 
L = 3.00; 
h = 250; %W/m^2K
T_in = 1150; %ºC
T_amb = 20;
ro = 7930;
c = 385;
k = 121;
v = 39.6*10^-6;
alpha = k/(ro*c);
L_x = 0.12;
L_y = 0.24;
L_z = 1.50;
Bi_x = h*L_x/k;
Bi_y = h*L_y/k;
Bi_z = h*L_z/k;

%% alinea b) - cálculo dos csi's
Bi = [Bi_x Bi_y Bi_z];
csi = zeros(2,50)

for l=1:3
    fun = @(csi)csi*tan(csi)-Bi(l);
    j=1;
    for i=1:50
        out = fzero(fun, i-1);
        if abs(fun(out))<0.05
            if j==1
                csi(l,j)=out;
                j = j+1;
            else if out ~= csi(l,j-1)
                    csi(l,j) = out;
                    j= j+1;
                end
            end
        end
    end
end
csi

%% Solução exata - 3D

ksi_x = [0.4782 3.2185 6.3224 9.4510 12.5861 15.7237 18.8627];
ksi_y = [0.6510 3.2911 6.3610 9.4771 12.6057 15.7395 18.8758];
ksi_z = [1.2011 3.8228 6.7156 9.7330 12.8039 15.9005 19.0112];

x = linspace(0, 1, 20);
y = linspace(0.5, 1, 20);
z = linspace(0, 1, 20);
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;