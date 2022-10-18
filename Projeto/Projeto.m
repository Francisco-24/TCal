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

%% alinea a) - lumped capacitance method
V_corpo = H*W*L;
A_corpo = H*W*2 + H*L*2 + W*L;
L_c = V_corpo/A_corpo;

Bi = h*L_c/k;

t_lcm = linspace(0, 8000, 200);
theta_star_lcm = @(t_lcm)exp(-h*A_corpo/(ro*V_corpo*c).*t_lcm);
figure()
plot(t_lcm, theta_star_lcm(t_lcm));

%% alinea b) - cálculo dos csi's
% Bi = [0.2479 0.4959 3.0992];
Bi = 3.0992
csi = 0


fun = @(csi)csi*tan(csi)-Bi;
j=1;
for i=1:50
    out = fzero(fun, i-1);
    if abs(fun(out))<0.05
        if j==1
            csi(j)=out;
            j = j+1;
        else if out ~= csi(j-1)
                csi(j) = out;
                j= j+1;
            end
        end
    end
end

csi
%% alinea b) - solução exata

alpha = k/(ro*c);
L_x = 0.12;
L_y = 0.24;
Bi_x = h*L_x/k
Bi_y = h*L_y/k
ksi_x = [0.4328 3.2039 6.3148 9.4459];
ksi_y = [0.6533 3.2923 6.3616 9.4775];

for i=1:4
    C_x(i) = 4*sin(ksi_x(i))/(2*ksi_x(i) + sin(2*ksi_x(i)));
    C_y(i) = 4*sin(ksi_y(i))/(2*ksi_y(i) + sin(2*ksi_y(i)));
end

t = linspace(0, 8000, 200);
theta_star_x = C_x(1)*exp(-ksi_x(1)^2*alpha.*t/(H/2)^2) + C_x(2)*exp(-ksi_x(2)^2*alpha.*t/(H/2)^2) + C_x(3)*exp(-ksi_x(3)^2*alpha.*t/(H/2)^2) +C_x(4)*exp(-ksi_x(4)^2*alpha.*t/(H/2)^2);
theta_star_y = C_y(1)*exp(-ksi_y(1)^2*alpha.*t/(H/2)^2)*cos(ksi_y(1)*0.5) + C_y(2)*exp(-ksi_y(2)^2*alpha.*t/(W)^2)*cos(ksi_y(2)*0.5) + C_y(3)*exp(-ksi_y(2)^2*alpha.*t/(H/2)^2)*cos(ksi_y(3)*0.5) + C_y(4)*exp(-ksi_y(4)^2*alpha.*t/(H/2)^2)*cos(ksi_y(4)*0.5);

figure()
plot(t, theta_star_x);
ylabel("$\theta*$ em x=0", 'Interpreter','latex', 'FontSize', 18)
xlabel("tempo [t]", 'FontSize', 12)
close

figure()
plot(t, theta_star_y);
ylabel("$\theta*$  em y=0.5", 'Interpreter','latex', 'FontSize', 18)
xlabel("tempo [t]", 'FontSize', 12)
close

figure()
plot(t, theta_star_x.*theta_star_y);
ylabel("$\theta*$  em 2D no centro do corpo", 'Interpreter','latex', 'FontSize', 18)
xlabel("tempo [t]", 'FontSize', 12)
close

x = [0 0.2 0.4 0.6 0.8 1];
y = [0.5 0.6 0.7 0.8 0.9 1];
theta_estrela_x = 0;
theta_estrela_y = 0;

theta_estrela = zeros(6,6);

for i=1:6
    for j=1:6
        for p=1:4
            theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha*500/(H/2)^2)*cos(ksi_x(p)*x(i));
            theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha*500/(W)^2)*cos(ksi_y(p)*y(j));
        end
        theta_estrela(i,j) = theta_estrela_x*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
    end
end

figure()
plot(y, theta_estrela, '.', 'markersize', 20)
legend({'x=0','x=0.2','x=0.4','x=0.6','x=0.8','x=1'},'Location','northeast','Orientation','vertical')
ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
xlabel("y [adimensional]", 'FontSize', 12)
close






