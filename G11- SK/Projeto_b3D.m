close all
clear all
clc
%% Dados
H = 0.24; 
W = 0.24; 
L = 3.00; 
h = 250;
T_in = 1150; 
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

for l=1:3
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

%% Solução exata - 3D

for i=1:20
    ksi_x(i)=csi(1,i+1);
    ksi_y(i)=csi(2,i+1);
    ksi_z(i)=csi(3,i+1);
end

x = [0 1];
y = [0 0.5 1];
z = linspace(0, 1, 40);
t = [0 250 500 1000 2000]
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;
theta_estrela = 0;

for i=1:2
    for j=1:3
        figure()
        for w=1:length(t)
            for k=1:length(z)
                for p=1:length(ksi_x)
                    C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                    C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                    C_z(p) = 4*sin(ksi_z(p))/(2*ksi_z(p) + sin(2*ksi_z(p)));
                    theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha*t(w)/(H/2)^2)*cos(ksi_x(p)*x(i));
                    theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha*t(w)/(W)^2)*cos(ksi_y(p)*y(j));
                    theta_estrela_z = theta_estrela_z + C_z(p)*exp(-ksi_z(p)^2*alpha*t(w)/(L/2)^2)*cos(ksi_z(p)*z(k));
                end
                theta_estrela(k) = theta_estrela_x*theta_estrela_y*theta_estrela_z;
                theta_estrela_x = 0;
                theta_estrela_y = 0;
                theta_estrela_z = 0;
            end
            plot(z, theta_estrela, 'LineWidth',1.5)
            legend(sprintf('t = %g s', t(1)),sprintf('t = %g s', t(2)),sprintf('t = %g s', t(3)),sprintf('t = %g s', t(4)),sprintf('t = %g s', t(5)),'Location','best','Orientation','vertical','FontSize', 10)
            ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
            xlabel("z*", 'FontSize', 20)
            title(['x* = ',num2str(x(i)),' y* =',num2str(y(j))])
            hold on
        end
    end
end

V_corpo = H*W*L;
A_corpo = H*W*2 + H*L*2 + W*L;
x = [0 1];
y = [0 0.5 1];
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;
theta_estrela = 0;
z = [0 1];
t = linspace(0, 6000, 300);
for k=1:length(z)
    for i=1:2
        for j=1:3
            for p=1:length(ksi_x)
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                C_z(p) = 4*sin(ksi_z(p))/(2*ksi_z(p) + sin(2*ksi_z(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(j));
                theta_estrela_z = theta_estrela_z + C_z(p)*exp(-ksi_z(p)^2*alpha.*t/(L/2)^2)*cos(ksi_z(p)*z(k));
            end
            theta_estrela_3D = theta_estrela_x.*theta_estrela_y.*theta_estrela_z; 
            theta_estrela_x = 0;
            theta_estrela_y = 0;
            theta_estrela_z = 0;
            figure()
            plot(t*alpha/(H/2)^2, theta_estrela_3D)
            title('Solução analítica 3D')
            legend(sprintf('x* = %g, y* = %g, z* = %g', x(i), y(j), z(k)),'Location','northeast','Orientation','vertical','FontSize', 15)
            ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
            hold on
        end
    end
end

