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

%% alinea a) - lumped capacitance method
V_corpo = H*W*L
A_corpo = H*W*2 + H*L*2 + W*L
L_c = V_corpo/A_corpo

Bi = h*L_c/k

t_lcm = linspace(0, 6000, 300);
theta_star_lcm = @(t_lcm)exp(-h*A_corpo/(ro*V_corpo*c).*t_lcm);
figure()
plot(t_lcm*alpha/(H/2)^2, theta_star_lcm(t_lcm),'LineWidth',2.0);
title('Método da Capacitância Global')
ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
xlabel("Fo", 'FontSize', 20)







