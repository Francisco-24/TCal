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
Bi_x = h*L_x/k
Bi_y = h*L_y/k

%% alinea b) - cálculo dos csi's
Bi = [Bi_x Bi_y];
csi = zeros(2,50)

for l=1:2
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

%% Solução exata - 2D

ksi_x = [0.4782 3.2185 6.3224  9.4510 12.5861 15.7237 18.8627];
ksi_y = [0.6510 3.2911 6.3610  9.4771 12.6057 15.7395 18.8758];

x = linspace(0, 1, 20);
y = linspace(0.5, 1, 20);
theta_estrela_x = 0;
theta_estrela_y = 0;

%theta_estrela = zeros(6,6);

for i=1:length(x)
    for j=1:length(y)
        for p=1:7
            C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
            C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
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

figure()
surf(x, y, theta_estrela);
shading interp
xlabel("x*", 'FontSize', 12)
ylabel("y*", 'FontSize', 12)
zlabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)

x=[0 1];
y=[0.5 1];
theta_estrela = 0;

t = linspace(0, 8000, 10);
for i=1:2
    for k=1:2
        for j=1:10
            for p=1:7
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha*t(j)/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha*t(j)/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela(j) = theta_estrela_x*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        end
        figure()
        plot(t, theta_estrela, '.', 'markersize', 20)
        legend(sprintf('x = %g, y = %g', x(i), y(k)),'Location','northeast','Orientation','vertical')
        ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
        xlabel("t", 'FontSize', 12)
        close
    end
end

t = linspace(0, 8000, 200);
for i=1:2
    for k=1:2
            for p=1:7
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela = theta_estrela_x.*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        figure()
        plot(t, theta_estrela)
        legend(sprintf('x = %g, y = %g', x(i), y(k)),'Location','northeast','Orientation','vertical')
        ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
        xlabel("t", 'FontSize', 12)
        hold on
    end
end


   


