close all
clear all
clc

%% pinos

D = 13
S_T = 28;
S_L = 17;
S_D = sqrt(S_L^2 + (S_T/2)^2);

%% 1º ensaio

V1 = 0.84495;
T_in1 = 26;
T_out1 = (34 + 33 +31)/3;

T_f1 = T_in1 + T_out1

if 2*(S_D - D)>= (S_T - D)
    V_max1 = S_T*V1/(S_T - D)
end



%% 2º ensaio

T_in2 = 26;
T_out2 = (32 + 31 + 29)/3;
V2 = 1.757;


T_f2 = T_in2 + T_out2

if 2*(S_D - D)>= (S_T - D)
    V_max2 = S_T*V2/(S_T - D)
end