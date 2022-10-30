close all
clear all
clc

% Lê ficheiro excel e cria variáveis com os valores obtidos para cada mesh

T_1 = readmatrix("Grid_independency.xlsx",'Sheet',1); % x = 0.0120 m; y = 0.2280 m
T_2 = readmatrix("Grid_independency.xlsx",'Sheet',2); % x = 0.0120 m; y = 0.1080 m
T_3 = readmatrix("Grid_independency.xlsx",'Sheet',3); % x = 0.0120 m; y = 0.0120 m
T_4 = readmatrix("Grid_independency.xlsx",'Sheet',4); % x = 0.0600 m; y = 0.2280 m
T_5 = readmatrix("Grid_independency.xlsx",'Sheet',5); % x = 0.0060 m; y = 0.1080 m
T_6 = readmatrix("Grid_independency.xlsx",'Sheet',6); % x = 0.0060 m; y = 0.0120 m
T_7 = readmatrix("Grid_independency.xlsx",'Sheet',7); % x = 0.1080 m; y = 0.2280 m
T_8 = readmatrix("Grid_independency.xlsx",'Sheet',8); % x = 0.1080 m; y = 0.0120 m
T_9 = readmatrix("Grid_independency.xlsx",'Sheet',9); % x = 0.1080 m; y = 0.0120 m

% Plots para cada ponto e mesh
figure()
plot(T_1(:,2),T_1(:,1),'DisplayName','12x12');
hold on
plot(T_1(:,5),T_1(:,4),'DisplayName','22x22');
hold on
plot(T_1(:,8),T_1(:,7),'DisplayName','42x42');
hold on
plot(T_1(:,11),T_1(:,10),'DisplayName','62x62');
hold on
plot(T_1(:,14),T_1(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.0120 m e y=0.2280 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_2(:,2),T_2(:,1),'DisplayName','12x12');
hold on
plot(T_2(:,5),T_2(:,4),'DisplayName','22x22');
hold on
plot(T_2(:,8),T_2(:,7),'DisplayName','42x42');
hold on
plot(T_2(:,11),T_2(:,10),'DisplayName','62x62');
hold on
plot(T_2(:,14),T_2(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.0120 m e y=0.1080 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_3(:,2),T_3(:,1),'DisplayName','12x12');
hold on
plot(T_3(:,5),T_3(:,4),'DisplayName','22x22');
hold on
plot(T_3(:,8),T_3(:,7),'DisplayName','42x42');
hold on
plot(T_3(:,11),T_3(:,10),'DisplayName','62x62');
hold on
plot(T_3(:,14),T_3(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.0120 m e y=0.0120 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_4(:,2),T_4(:,1),'DisplayName','12x12');
hold on
plot(T_4(:,5),T_4(:,4),'DisplayName','22x22');
hold on
plot(T_4(:,8),T_4(:,7),'DisplayName','42x42');
hold on
plot(T_4(:,11),T_4(:,10),'DisplayName','62x62');
hold on
plot(T_4(:,14),T_4(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.0600 m e y=0.2280 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_5(:,2),T_5(:,1),'DisplayName','12x12');
hold on
plot(T_5(:,5),T_5(:,4),'DisplayName','22x22');
hold on
plot(T_5(:,8),T_5(:,7),'DisplayName','42x42');
hold on
plot(T_5(:,11),T_5(:,10),'DisplayName','62x62');
hold on
plot(T_5(:,14),T_5(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.0600 m e y=0.1080 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_6(:,2),T_6(:,1),'DisplayName','12x12');
hold on
plot(T_6(:,5),T_6(:,4),'DisplayName','22x22');
hold on
plot(T_6(:,8),T_6(:,7),'DisplayName','42x42');
hold on
plot(T_6(:,11),T_6(:,10),'DisplayName','62x62');
hold on
plot(T_6(:,14),T_6(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.0600 m e y=0.0120 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_7(:,2),T_7(:,1),'DisplayName','12x12');
hold on
plot(T_7(:,5),T_7(:,4),'DisplayName','22x22');
hold on
plot(T_7(:,8),T_7(:,7),'DisplayName','42x42');
hold on
plot(T_7(:,11),T_7(:,10),'DisplayName','62x62');
hold on
plot(T_7(:,14),T_7(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.1080 m e y=0.2280 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_8(:,2),T_8(:,1),'DisplayName','12x12');
hold on
plot(T_8(:,5),T_8(:,4),'DisplayName','22x22');
hold on
plot(T_8(:,8),T_8(:,7),'DisplayName','42x42');
hold on
plot(T_8(:,11),T_8(:,10),'DisplayName','62x62');
hold on
plot(T_8(:,14),T_8(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.1080m e y=0.1080m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')

figure()
plot(T_9(:,2),T_9(:,1),'DisplayName','12x12');
hold on
plot(T_9(:,5),T_9(:,4),'DisplayName','22x22');
hold on
plot(T_9(:,8),T_9(:,7),'DisplayName','42x42');
hold on
plot(T_9(:,11),T_9(:,10),'DisplayName','62x62')
hold on
plot(T_9(:,14),T_9(:,13),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.1080 m e y=0.0120 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')


% Temperatura obtida para a solução analítica 2D para as coordenadas
% escolhidas

T = [969.121179933070	889.031901896929	607.716478082218	389.813461374913	167.202223429503	78.5975530176150	43.3262334287646	29.2855953530158;
1052.59913471148	994.989580820830	708.525756700838	454.502223597047	192.965583173422	88.8533117695969	47.4087967850567	30.9107626328236;
1055.55081156739	1008.15696304844	737.331590807199	473.569757660760	200.566075252144	91.8788824216115	48.6132014693804	31.3902062847131;
1025.91837939899	947.537925154420	648.304038514685	415.352707329495	177.367980444546	82.6442886699369	44.9371385875890	29.9268567680921;
1114.39181225628	1060.62901588277	756.075181979154	484.508862930838	204.910552818181	93.6083081332625	49.3016429778660	31.6642577852065;
1117.52012301594	1074.68286864059	786.870340149415	504.893197200769	213.035933406511	96.8428241080015	50.5892236162165	32.1768117227558;
1045.60485543228	972.120155172193	665.967176759620	426.467039113580	181.791979355970	84.4053728758843	45.6381825588840	30.2059249966595;
1135.80977083447	1088.20846157211	776.768026368152	497.567343430779	210.108840818494	95.6776177484295	50.1253838447871	31.9921686066998;
1138.99930466755	1102.63477892869	808.428910528845	518.524731204767	218.462645732892	99.0030638800581	51.4491615247870	32.5191316893329];

% Calcula erro quadrático médio para cada mesh e instante
% Erro quadrático médio para t = 60s

for m=1:3:13
    EQM60(m) = ((T_1(1,m)-T(1,1))^2+(T_2(1,m)-T(2,1))^2+(T_3(1,m)-T(3,1))^2+(T_4(1,m)-T(4,1))^2+(T_5(1,m)-T(5,1))^2+(T_6(1,m)-T(6,1))^2+(T_7(1,m)-T(7,1))^2+(T_8(1,m)-T(8,1))^2+(T_9(1,m)-T(9,1))^2)/9;
end 

% Erro quadrático médio para t = 120s
for m=1:3:13
    EQM120(m) = ((T_1(2,m)-T(1,2))^2+(T_2(2,m)-T(2,2))^2+(T_3(2,m)-T(3,2))^2+(T_4(2,m)-T(4,2))^2+(T_5(2,m)-T(5,2))^2+(T_6(2,m)-T(6,2))^2+(T_7(2,m)-T(7,2))^2+(T_8(2,m)-T(8,2))^2+(T_9(2,m)-T(9,2))^2)/9;
end

% Erro quadrático médio para t = 500s
for m=1:3:13
    EQM500(m) = ((T_1(3,m)-T(1,3))^2+(T_2(3,m)-T(2,3))^2+(T_3(3,m)-T(3,3))^2+(T_4(3,m)-T(4,3))^2+(T_5(3,m)-T(5,3))^2+(T_6(3,m)-T(6,3))^2+(T_7(3,m)-T(7,3))^2+(T_8(3,m)-T(8,3))^2+(T_9(3,m)-T(9,3))^2)/9;
end

% Erro quadrático médio para t = 1000s
for m=1:3:13
    EQM1000(m) = ((T_1(4,m)-T(1,4))^2+(T_2(4,m)-T(2,4))^2+(T_3(4,m)-T(3,4))^2+(T_4(4,m)-T(4,4))^2+(T_5(4,m)-T(5,4))^2+(T_6(4,m)-T(6,4))^2+(T_7(4,m)-T(7,4))^2+(T_8(4,m)-T(8,4))^2+(T_9(4,m)-T(9,4))^2)/9;
end

%Erro quadrático médio para t = 2000s
for m=1:3:13
    EQM2000(m) = ((T_1(5,m)-T(1,5))^2+(T_2(5,m)-T(2,5))^2+(T_3(5,m)-T(3,5))^2+(T_4(5,m)-T(4,5))^2+(T_5(5,m)-T(5,5))^2+(T_6(5,m)-T(6,5))^2+(T_7(5,m)-T(7,5))^2+(T_8(5,m)-T(8,5))^2+(T_9(5,m)-T(9,5))^2)/9;
end

%Erro quadrático médio para t = 3000s
for m=1:3:13
    EQM3000(m) = ((T_1(6,m)-T(1,6))^2+(T_2(6,m)-T(2,6))^2+(T_3(6,m)-T(3,6))^2+(T_4(6,m)-T(4,6))^2+(T_5(6,m)-T(5,6))^2+(T_6(6,m)-T(6,6))^2+(T_7(6,m)-T(7,6))^2+(T_8(6,m)-T(8,6))^2+(T_9(6,m)-T(9,6))^2)/9;
end

%Erro quadrático médio para t = 3000s
for m=1:3:13
    EQM4000(m) = ((T_1(7,m)-T(1,7))^2+(T_2(7,m)-T(2,7))^2+(T_3(7,m)-T(3,7))^2+(T_4(7,m)-T(4,7))^2+(T_5(7,m)-T(5,7))^2+(T_6(7,m)-T(6,7))^2+(T_7(7,m)-T(7,7))^2+(T_8(7,m)-T(8,7))^2+(T_9(7,m)-T(9,7))^2)/9;
end

%Erro quadrático médio para t = 3000s
for m=1:3:13
    EQM5000(m) = ((T_1(8,m)-T(1,8))^2+(T_2(8,m)-T(2,8))^2+(T_3(8,m)-T(3,8))^2+(T_4(8,m)-T(4,8))^2+(T_5(8,m)-T(5,8))^2+(T_6(8,m)-T(6,8))^2+(T_7(8,m)-T(7,8))^2+(T_8(8,m)-T(8,8))^2+(T_9(8,m)-T(9,8))^2)/9;
end

EQM60 = nonzeros(EQM60);
EQM120 = nonzeros(EQM120);
EQM500 = nonzeros(EQM500);
EQM1000 = nonzeros(EQM1000);
EQM2000 = nonzeros(EQM2000);
EQM3000 = nonzeros(EQM3000);
EQM4000 = nonzeros(EQM4000);
EQM5000 = nonzeros(EQM5000);

x = [12 22 42 62 82];

figure()
plot(x,EQM60)
title(['Erros quadráticos médios em função da malha para TIME = 60s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);


figure()
plot(x,EQM120)
title(['Erros quadráticos médios em função da malha para TIME = 120s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);

figure()
plot(x,EQM500)
title(['Erros quadráticos médios em função da malha para TIME = 500s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);

figure()
plot(x,EQM1000)
title(['Erros quadráticos médios em função da malha para TIME = 1000s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);

figure()
plot(x,EQM2000)
title(['Erros quadráticos médios em função da malha para TIME = 2000s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);

figure()
plot(x,EQM3000)
title(['Erros quadráticos médios em função da malha para TIME = 3000s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);

figure()
plot(x,EQM4000)
title(['Erros quadráticos médios em função da malha para TIME = 4000s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);

figure()
plot(x,EQM5000)
title(['Erros quadráticos médios em função da malha para TIME = 5000s'],'FontSize',24);
xlabel('Número de nós da malha segundo x e y','Interpreter','latex','FontSize',24);
ylabel('Erro quadrático médio','Interpreter','latex','FontSize',24);
