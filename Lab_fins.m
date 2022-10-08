close all
clear all
clc

%% pinos

L = 0.067;
A = 0.120*0.070;
D = 13;
S_T = 28;
S_L = 17;
S_D = sqrt(S_L^2 + (S_T/2)^2);

k_al = 237; 
L_c = (L + D/4*10^-3); 
A_f = pi*D*10^-3*L;
A_t = 17*A_f + 0.110*0.100 - 17*pi*(6.5*10^-3)^2;
A_conv = A_t;


%% 1º ensaio

V1 = 0.84495;
T_in1 = 23.5;
T_out1 = (31.5 + 30.5 + 28.5)/3;
T_s1 = 49.5;
T_f1 = (T_in1 + T_out1)/2 + 273.15;

if 2*(S_D - D)>= (S_T - D)
    V_max1 = S_T*V1/(S_T - D);
end

%propriedades para a temperatura de filme térmico
ro1 =  1.1614 + ((T_f1-300)/(350-300)*(0.9950-1.1614));
cp1 = 1007 + ((T_f1-300)/(350-300)*(1009-1007));
niu1 = (184.6 + ((T_f1-300)/(350-300)*(208.2-184.6)))*10^-7;
v1 = (15.89 + ((T_f1-300)/(350-300)*(20.92-15.89)))*10^-6;
k1 = (26.3 + ((T_f1-300)/(350-300)*(30-26.3)))*10^-3;
alpha1 = (22.5 + ((T_f1-300)/(350-300)*(29.9-22.5)))*10^-6;
Pr1 = 0.707 + ((T_f1-300)/(350-300)*(0.700-0.707));
%numero de Prandl para a temperatura da base
Prs1 = 0.707 + ((T_s1+275.15-300)/(350-300)*(0.700-0.707));

Re_max1 = V_max1*D*10^-3/v1;
if S_T/S_L<2
    C1 = 0.35*(S_T/S_L)^(1/5);
    m = 0.60;
end

Nu1 = 0.92*C1*(Re_max1^m)*Pr1^0.36*(Pr1/Prs1)^0.25;

h_teorico1 = Nu1*k1/D*10^3

mf1 = ro1 * A * V1;
q1 = mf1*cp1*(T_out1-T_in1);
DeltaT1 = ((T_s1 - T_in1) - (T_s1 - T_out1))/log((T_s1 - T_in1)/(T_s1 - T_out1));

h_exp1 = q1/(A_conv*DeltaT1)

m1 = (h_teorico1*pi*(D*10^-3)/(k_al*pi*(D*10^-3)^2/4))^0.5;
eta_teo1 = tanh(m1*L_c)/(m1*L_c)
eta_teo_total = 1 - 17*A_f/A_t*(1-eta_teo1)

%eta_exp = q1/(h_exp1*A_f*(T_s1-T_in1))
eta_exp_total = q1/(h_exp1*A_t*(T_s1-T_in1))

x1 = [10, 36, 62];
y1 = [0.731, 0.615, 0.538];

y1_x = @(x)(cosh(m1.*(L-x)) + h_teorico1./(m1.*k_al).*sinh(m1.*(L-x)))/(cosh(m1.*L) + h_teorico1./(m1.*k_al).*sinh(m1*L));
x = linspace(0, 0.067, 100);
figure()
plot(x*1000, y1_x(x))
hold on
plot(x1, y1, '.', 'markersize', 10)
ylim([0 1])
ylabel("$\frac{\theta}{\theta_b}$", 'Interpreter','latex')
xlabel("distance [mm]")

%y1_x_exp = @(x)(cosh(m1.*(L-x)) + h_exp1./(m1.*k_al).*sinh(m1.*(L-x)))/(cosh(m1.*L) + h_exp1./(m1.*k_al).*sinh(m1*L));
%x = linspace(0, 0.067, 100);
%figure()
%plot(x*1000, y1_x_exp(x))
%hold on
%plot(x1, y1, '.', 'markersize', 10)
%ylim([0 1])
%ylabel("$\frac{\theta}{\theta_b}$", 'Interpreter','latex')
%xlabel("distance [mm]")


%% 2º ensaio

T_in2 = 23.5;
T_out2 = (29.5 + 28.5 + 26-5)/3;
V2 = 1.757;
T_f2 = (T_in2 + T_out2)/2 + 273.15;
T_s2 = 39.5;

if 2*(S_D - D)>= (S_T - D)
    V_max2 = S_T*V2/(S_T - D);
end

ro2 =  1.1614 + ((T_f2-300)/(350-300)*(0.9950-1.1614));
cp2 = 1007 + ((T_f2-300)/(350-300)*(1009-1007));
niu2 = (184.6 + ((T_f2-300)/(350-300)*(208.2-184.6)))*10^-7;
v2 = (15.89 + ((T_f2-300)/(350-300)*(20.92-15.89)))*10^-6;
k2 = (26.3 + ((T_f2-300)/(350-300)*(30-26.3)))*10^-3;
alpha2 = (22.5 + ((T_f2-300)/(350-300)*(29.9-22.5)))*10^-6;
Pr2 = 0.707 + ((T_f2-300)/(350-300)*(0.700-0.707));
Prs2 = 0.707 + ((T_s2+275.15-300)/(350-300)*(0.700-0.707));

Re_max2 = V_max2*D*10^-3/v2;
if S_T/S_L<2
    C2 = 0.35*(S_T/S_L)^(1/5);
    m = 0.60;
end

Nu2 = 0.92*C2*(Re_max2^m)*Pr2^0.36*(Pr2/Prs2)^0.25;

h_teorico2 = Nu2*k2/D*10^3

mf2 = ro2 * A * V2;
q2 = mf2*cp2*(T_out2-T_in2);
DeltaT2 = ((T_s2 - T_in2) - (T_s2 - T_out2))/log((T_s2 - T_in2)/(T_s2 - T_out2));

h_exp2 = q2/(A_conv*DeltaT2)

m2 = sqrt(h_teorico2*pi*(D*10^-3)/(k_al*pi*(D*10^-3)^2/4));
eta_teo2 = tanh(m2*L_c)/(m2*L_c)
eta_teo_total2 = 1 - 17*A_f/A_t*(1-eta_teo2)

%eta_exp = q1/(h_exp1*A_f*(T_s1-T_in1))
eta_exp_total2 = q2/(h_exp2*A_t*(T_s2-T_in2))

x2 = [10, 36, 62];
y2 = [0.75, 0.625, 0.562];

y2_x = @(x)(cosh(m2.*(L-x)) + h_teorico2./(m2.*k_al).*sinh(m2.*(L-x)))/(cosh(m2.*L) + h_teorico2./(m2.*k_al).*sinh(m2*L));
x = linspace(0, 0.067, 100);
figure()
plot(x*1000, y2_x(x))
hold on
plot(x2, y2, '.', 'markersize', 8)
ylim([0 1])
ylabel("$\frac{\theta}{\theta_b}$", 'Interpreter','latex')
xlabel("distance [mm]")

%y2_x_exp = @(x)(cosh(m2.*(L-x)) + h_exp2./(m2.*k_al).*sinh(m2.*(L-x)))/(cosh(m2.*L) + h_exp2./(m2.*k_al).*sinh(m2*L));
%x = linspace(0, 0.067, 100);
%figure()
%plot(x*1000, y2_x_exp(x))
%hold on
%plot(x2, y2, '.', 'markersize', 8)
%ylim([0 1])
%ylabel("$\frac{\theta}{\theta_b}$", 'Interpreter','latex')
%xlabel("distance [mm]")

