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
Bi = [Bi_x Bi_y Bi_z]


for l=1:3
    fun = @(csi)csi*tan(csi)-Bi(l);
    j=1;
    for i=1:200
        out = fzero(fun, i-1);
        if abs(fun(out))<0.05
            if j==1
                csi(l,j)=out;
                j = j+1;
            else if out ~= csi(l,j-1) && out-csi(l,j-1)>0.5
                    csi(l,j) = out;
                    j= j+1;
                end
            end
        end
    end
end
csi;


%% Solução exata - 3D

ksi_x = [0.4782 3.2185 6.3224 9.4510 12.5861 15.7237 18.8627 22.0024 25.1426 28.2831]
ksi_y = [0.6510 3.2911 6.3610 9.4771 12.6057 15.7395 18.8758 22.0137 25.1525 28.2919]
ksi_z = [1.2011 3.8228 6.7156 9.7330 12.8039 15.9005 19.0112 22.1303 25.2548 28.3831]

for i=1:20
    ksi_x(i)=csi(1,i+1);
    ksi_y(i)=csi(2,i+1);
    ksi_z(i)=csi(3,i+1);
end

ksi_x 
ksi_y 
ksi_z 

x = linspace(0, 1, 6);
y = linspace(0, 1, 20);
z = linspace(0, 1, 5);
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;

%for em que para cada gráfico existe um único z
for k=1:length(z)
%     figure()
    for i=1:length(x)
        for j=1:length(y)
            for p=1:length(ksi_x)
                C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
                C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
                C_z(p) = 4*sin(ksi_z(p))/(2*ksi_z(p) + sin(2*ksi_z(p)));
                theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha*500/(H/2)^2)*cos(ksi_x(p)*x(i));
                theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha*500/(W)^2)*cos(ksi_y(p)*y(j));
                theta_estrela_z = theta_estrela_z + C_z(p)*exp(-ksi_z(p)^2*alpha*500/(L/2)^2)*cos(ksi_z(p)*z(k));
            end
            theta_estrela(j,i) = theta_estrela_x*theta_estrela_y*theta_estrela_z;
            theta_teste(j) = theta_estrela_x*theta_estrela_y;
            theta_estrela_x = 0;
            theta_estrela_y = 0;
            theta_estrela_z = 0;
        end
%         plot(y, theta_teste, '.', 'markersize', 20)
%         legend({'x=0','x=0.2','x=0.4','x=0.6','x=0.8','x=1'},'Location','northeast','Orientation','vertical')
%         ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
%         xlabel("y [adimensional]", 'FontSize', 12)
%         hold on
    end

figure()
plot(y, theta_estrela, '.', 'markersize', 20)
legend({'x=0','x=0.2','x=0.4','x=0.6','x=0.8','x=1'},'Location','northeast','Orientation','vertical')
ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
xlabel("y [adimensional]", 'FontSize', 12)
title( ['z=', num2str(z(k))])
end

x = [0 1];
y = [0 0.5 1];
z = linspace(0, 1, 40);
t = linspace(1, 3000, 5);
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;
theta_estrela = 0;

%for para os pontos x e y fixos e variar o z e varia o t

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
            plot(z, theta_estrela)
            %legend({'t=1 s','t=250 s','t=500 s','t=750 s','t=1000 s'},'Location','best','Orientation','vertical')
            ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
            xlabel("z [adimensional]", 'FontSize', 12)
            title(['x = ',num2str(x(i)),' y =',num2str(y(j))])
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
t = linspace(0, 4000, 200);
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
            theta_estrela_2D = theta_estrela_x.*theta_estrela_y;
            theta_estrela_3D = theta_estrela_x.*theta_estrela_y.*theta_estrela_z;
            theta_star_lcm = exp(-h*A_corpo/(ro*V_corpo*c).*t);
            erro_2D_3_D = abs(theta_estrela_3D - theta_estrela_2D);
            erro_LCM_3D = abs(theta_star_lcm - theta_estrela_3D);
            theta_estrela_x = 0;
            theta_estrela_y = 0;
            theta_estrela_z = 0;
            figure()
            plot(t*alpha/(H/2)^2, theta_estrela_2D,'--',t*alpha/(H/2)^2, theta_estrela_3D, '+', t*alpha/(H/2)^2, theta_star_lcm)
            legend(sprintf('x = %g, y = %g 2D', x(i), y(j)),sprintf('x = %g, y = %g, z = %g 3D', x(i), y(j), z(k)), sprintf('LCM'),'Location','northeast','Orientation','vertical')
            ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
            xlabel("Fo", 'FontSize', 12)
            hold on
        end
    end
end

%% erros
V_corpo = H*W*L;
A_corpo = H*W*2 + H*L*2 + W*L;
x = [0 1];
y = [0 0.5 1];
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;
theta_estrela = 0;
z = [0 1];
t = linspace(0, 4000, 200);
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
            theta_estrela_2D = theta_estrela_x.*theta_estrela_y;
            theta_estrela_3D = theta_estrela_x.*theta_estrela_y.*theta_estrela_z;
            theta_star_lcm = exp(-h*A_corpo/(ro*V_corpo*c).*t);
            erro_2D_3_D = abs(theta_estrela_3D - theta_estrela_2D);
            erro_LCM_3D = abs(theta_star_lcm - theta_estrela_3D);
            theta_estrela_x = 0;
            theta_estrela_y = 0;
            theta_estrela_z = 0;
            figure()
            plot(t,erro_2D_3_D,'--', t,erro_LCM_3D,'+')
            title(sprintf('x = %g, y = %g, z = %g 3D', x(i), y(j), z(k)))
            legend('ERRO 2D','ERRO LCM','Location','northeast','Orientation','vertical')
            ylabel("ERRO", 'Interpreter','latex', 'FontSize', 18)
            xlabel("t", 'FontSize', 12)
            
        end
    end
end

