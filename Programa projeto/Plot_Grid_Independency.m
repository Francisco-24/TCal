close all
clear all
clc


T_1 = readmatrix("Grid_independency.xlsx",'Sheet',1);
T_2 = readmatrix("Grid_independency.xlsx",'Sheet',2);
T_3 = readmatrix("Grid_independency.xlsx",'Sheet',3);
T_4 = readmatrix("Grid_independency.xlsx",'Sheet',4);
T_5 = readmatrix("Grid_independency.xlsx",'Sheet',5);
T_6 = readmatrix("Grid_independency.xlsx",'Sheet',6);
T_7 = readmatrix("Grid_independency.xlsx",'Sheet',7);
T_8 = readmatrix("Grid_independency.xlsx",'Sheet',8);
T_9 = readmatrix("Grid_independency.xlsx",'Sheet',9);

figure()
plot(T_1(:,2),T_1(:,1),'DisplayName','12x12');
hold on
plot(T_1(:,5),T_1(:,4),'DisplayName','22x22');
hold on
plot(T_1(:,8),T_1(:,7),'DisplayName','42x42');
hold on
plot(T_1(:,11),T_1(:,10),'DisplayName','82x82');
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
plot(T_2(:,11),T_2(:,10),'DisplayName','82x82');
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
plot(T_3(:,11),T_3(:,10),'DisplayName','82x82');
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
plot(T_4(:,11),T_4(:,10),'DisplayName','82x82');
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
plot(T_5(:,11),T_5(:,10),'DisplayName','82x82');
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
plot(T_6(:,11),T_6(:,10),'DisplayName','82x82');
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
plot(T_7(:,11),T_7(:,10),'DisplayName','82x82');
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
plot(T_8(:,11),T_8(:,10),'DisplayName','82x82');
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
plot(T_9(:,11),T_9(:,10),'DisplayName','82x82');
title('Temperatura em função tempo para x=0.1080 m e y=0.0120 m','FontSize',14);
xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
legend ('Location','best')