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
Bi_x = h*L_x/k;
Bi_y = h*L_y/k;

%% alinea b) - cálculo dos csi's
Bi = [Bi_x Bi_y];

for l=1:2
    fun = @(csi)csi*tan(csi)-Bi(l);
    j=1;
    for i=1:200
        out = fzero(fun, i-1);
        if abs(fun(out))<0.05
            if j==1
                csi(l,j)=out;
                j = j+1;
            elseif out ~= csi(l,j-1) && out-csi(l,j-1)>0.5
                    csi(l,j) = out;
                    j= j+1;
            end
        end
    end
end
csi;

%% Solução exata - 2D

for i=1:20
    ksi_x(i)=csi(1,i+1);
    ksi_y(i)=csi(2,i+1);
end

theta_estrela_x = 0;
theta_estrela_y = 0;
x=[0 1];
y=[0 0.5 1];
theta_estrela = 0;
t = linspace(0, 6000, 600);
for i=1:2
    for k=1:3
            for p=1:length(ksi_x)
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela = theta_estrela_x.*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        figure()
        plot(t*alpha/(H/2)^2, theta_estrela,'LineWidth',1.5)
        title('Solução analítica 2D')
        legend(sprintf('x* = %g, y* = %g', x(i), y(k)),'Location','northeast','Orientation','vertical','FontSize', 15)
        ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
        xlabel("Fo", 'FontSize', 20)
        hold on
    end
end

t = linspace(0, 6000, 600);
for i=1:2
    for k=1:2
            for p=1:length(ksi_x)
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela_20 = theta_estrela_x.*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        for p=1:10
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela_10 = theta_estrela_x.*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        for p=1:5
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela_5 = theta_estrela_x.*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        for p=1:1
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(k));
            end
        theta_estrela_1 = theta_estrela_x.*theta_estrela_y;
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        figure()
        plot(t*alpha/(H/2)^2, theta_estrela_20,t*alpha/(H/2)^2, theta_estrela_10,t*alpha/(H/2)^2,theta_estrela_5,t*alpha/(H/2)^2,theta_estrela_1,'LineWidth',1.5)
        title('Solução analítica 2D')
        legend(sprintf('x* = %g, y* = %g 20 termos', x(i), y(k)),sprintf('x* = %g, y* = %g 10 termos', x(i), y(k)),sprintf('x* = %g, y* = %g 5 termos', x(i), y(k)),sprintf('x* = %g, y* = %g 1 termo', x(i), y(k)),'Location','northeast','Orientation','vertical','FontSize', 15)
        ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
        xlabel("Fo", 'FontSize', 20)
        hold on
    end
end

%% Estudo da independência da malha

x_ad = [0.108 0.060 0.0120]/L_x; 
y_ad = [0.2280 0.1080 0.0120]/L_y; 
t = [60 120 500 1000 2000 3000 4000 5000];
h=1;
theta_estrela_x = 0;
theta_estrela_y = 0;

for i=1:3
    for k=1:3
        for j=1:8
            for p=1:length(ksi_x)
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha*t(j)/(H/2)^2)*cos(ksi_x(p)*x_ad(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha*t(j)/(W)^2)*cos(ksi_y(p)*y_ad(k));
            end
        T(h,j) = theta_estrela_x*theta_estrela_y*(T_in-T_amb)+T_amb;    
        theta_estrela_x = 0;
        theta_estrela_y = 0;      
        end
        h=h+1;
    end
end




   


