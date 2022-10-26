%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                   Programa TEACH C (vers�o MATLAB)
%                               CONTROL
%
%     Este programa permite resolver problemas de condu��o de calor 
%               bidimensional com as seguintes variantes:
%
%    - Coordenadas cartesianas ou coordenadas cil�ndricas (axissim�trico)
%    - Regime estacion�rio ou transiente (n�o estacion�rio)
%    - Condutividade t�rmica (k) uniforme ou vari�vel 
%      (com a temperatura ou com o espa�o)
%    - Com a utiliza��o dos tr�s tipos de condi��es de fronteira existentes
%
%          - Temperatura imposta
%          - Fluxo imposto
%          - Convec��o na fronteira do corpo
%
%                                 FS
%                              09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% Apaga vari�veis e limpa a consola
clear all
close all
clc

% Utiliza o formato longo (15 d�gitos) para maior precis�o
format long

% Cria ficheiro para escrever os resultados
fid=fopen ('RESULTS.txt','w');

%---------------------- Capitulo 0 - Preliminares ------------------------%

% Define o n�mero de n�s da malha segundo x e y
%%%%% Alter�vel %%%%%
IT=12;
JT=12;

% Define o n�mero de pontos da malha segundo x e y
%%%%% Alter�vel %%%%%
NI=12;
NJ=12;

% Constante
GREAT=1.0E30;


%------------ Capitulo 1 - Par�metros e �ndices de controlo --------------%

% Define o tipo de coordenadas (cartesianas- 0 0; 
% cil�ndricas vertical 1 0; cil�ndricas horizontais 0 1)
%%%%% Alter�vel %%%%%
INCYLX=0;
INCYLY=0;

% Vari�vel interm�dia
NIM1=NI-1;
NJM1=NJ-1;

% Dimens�es totais do dom�nio de solu��o, largura e altura [m]
%%%%% Alter�vel %%%%%
W=0.24;
H=0.24;


% C�lculo as abcissas na dire��o xx
DX(1)=W/(NIM1-1);
X(1)=-0.5*DX(1);
for I=2:NI
    DX(I)=DX(1);
    X(I)=X(I-1)+DX(I-1);
end

% Calcula as abcissas na dire��o yy
DY(1)=H/(NJM1-1);
Y(1)=-0.5*DY(1);
for J=2:NJ
    DY(J)=DY(1);
    Y(J)=Y(J-1)+DY(J-1);
end

% Estabelece limites do Dom�nio 
for I=1:NI
    JS (I)=2;
    JN (I)=NJ-1;   
end

% Estabelece ponto monitor segundo x e y (IMON e JMON)
%%%%% Alter�vel %%%%%
IMON=2;
JMON=11;

% Propriedades do material (Tcond=cond,CV=calor esp,DENSIT=dens), 
% Meio homog�neo
% A�o %%%%% Alter�vel %%%%%
for I=1:NI
    for J=1:NJ
        TCON (I,J)=121;
        CV (I,J)=385;
        DENSIT (I,J)=7930;
        if I==1 && J==1
           BK=TCON (I,J);
            
        end
    end
end

% Temperatura inicial
TINIC=1150;

% Par�metros de controlo do programa
%%%%% Alter�vel %%%%%
% N�mero m�ximo de itera��es
MAXIT=300;
% N�mero m�ximo de intera��es no tempo 
MAXSTP=300;
% O output dever� conter os valores de T em intervalos de
NITPRI=300;
% "NITPRI" para "NSTPRI" itera��es no tempo
NSTPRI=1;    

% Factor de sub relaxa��o, M�ximo res�duo e intervalo de tempo (s)
%%%%% Alter�vel %%%%%
URFT=1;
SORMAX=0.001;
DT=20;

% Selecciona o Regime---Estacion�rio->INTIME=0, Transiente->INTIME=1
%%%%% Alter�vel %%%%%
INTIME=1;

if INTIME==0
    MAXSTP=1;
end

% Indica se as propriedades s�o constantes --- constantes->INPRO=0,
% vari�veis->INPRO=1
%%%%% Alter�vel %%%%%
INPRO=0;


%------------------ Capitulo 2 - Opera\[CCedilla]\[OTilde]es iniciais ----------------------%

% Calcula dimens�es da malha e anula vectores/matrizes
% Chama fun��o INIT
[RX,DXEP,DXPW,SEW,XU,RU,RY,DYPS,DYNP,SNS,YV,RV,AN,AS,AE,AW,SU,SP,GAMH,TOLD,T,X,Y]=INIT (INCYLX,INCYLY,NI,NJ,NIM1,NJM1,X,Y,TINIC);
TIME=0.0;
 
% Imp�e valores de fronteira e inicializa vari�vel dependente

% Valores de fronteira
%%%%% Alter�vel %%%%%
TTOP=20;
TBOT=1150;
TLEFT=0;
TRIGHT=0;

% Uma vez que de acordo com o enunciado n�o h� temperatura imposta em
% nenhuma das faces estes valores n�o ir�o ter import�ncia nos resultados,
% estando apenas a adicionar tempo de computa��o

% for I=2:NIM1
%     T (I,1)=TBOT;
%     T (I,NJ)=TTOP;
% end
%  
% for J=2:NJM1
%     T (1,J)=TLEFT;
%     T (NI,J)=TRIGHT;
% end
 

% Inicializa vari�vel dependente
% Inicializa campo de propriedades do material
% Chama fun��o PROPS
[GAMH]=PROPS(NI,NJ,TCON,GAMH);
 
 
% C�lculo do factor de normaliza��o do res�duo
AK=BK;
SNORM=AK*(TTOP-TBOT)*W/H;
SNORM=abs(SNORM);

% Escreve as especifica��es do problema
fprintf(fid,'CONDUCTION IN RECTANGULAR BAR WITH PRESCRIBED SURFACE TEMPERATURE \r \n \r \n');
fprintf(fid,'HEIGHT, H [M]--------------------------------= %10.3f \r \n',H);
fprintf(fid,'WEIGHT, W [M]--------------------------------= %10.3f \r \n',W);
fprintf(fid,'SPECIFIC HEAT, CV [J/KG.K]-------------------= %10.3f \r \n',CV(1,1));
fprintf(fid,'THERMAL CONDUCTIVITY, TCON [W/M.K]-----------= %10.3f \r \n',TCON(1,1));
fprintf(fid,'DENSIT, DENSIT [KG/M3] ----------------------= %10.2f \r \n',DENSIT(1,1));
fprintf(fid,'INITIAL TIME STEP, DT [S] -------------------= %10.1f \r \n',DT);
fprintf(fid,'SOURCE NORMALIZATION FACTOR, SNORM ----------= %8.3E \r \n',SNORM);
fprintf(fid,'NUMBER OF NODES IN X DIRECTION, NI ----------= %10d \r \n',NI);
fprintf(fid,'NUMBER OF NODES IN Y DIRECTION, NJ ----------= %10d \r \n',NJ);
fprintf(fid,'\r \n');

% Chama a fun��o PRINT para imprimir o campo de temperaturas inicial
PRINT(1,1,NI,NJ,X,Y,T,fid)


%-------------- Capitulo 3 - Itera��es no tempo e no espa�o --------------%

% Indica os pontos de controlo
% Imprime o r�tulo das informa��es das itera��es no ponto monitor
fprintf (fid,'NITER        SOURCE        T (%d,%d)    TIME(s)     DT(s)    NSTEP',IMON,JMON);
fprintf (fid,'\r \n');

Temperatura_monitor = zeros(MAXSTP,2);
Temperatura_centrogeo = zeros(MAXSTP,2);
Temperatura_centrosurfN = zeros(MAXSTP,2);
Temperatura_centrosurfO = zeros(MAXSTP,2);
Tmesh = zeros(6,2,9);
c = 1;
Pts_grid = [2,11;2,6;2,2;4,11;4,6;4,2;6,11;6,6;6,2];

% Itera��es no tempo
for NSTEP=1:MAXSTP
    TIME=TIME+DT;
    for I=1:NI
        for J=1:NJ
            TOLD(I,J)=T(I,J);
        end
    end
  
    % Itera��es no espa�o
    for NITER=1:MAXIT
        
        % Chama a fun��o CALCT para o c�lculo das temperaturas
        [AN,AS,AE,AW,SU,SP,GAMH,CV,DENSIT,TOLD,T,RESORT]=CALCT(NI,NJ,NIM1,NJM1,RX,DXEP,SEW,RU,RY,DYNP,SNS,RV,AN,AS,AE,AW,SU,SP,GAMH,CV,DENSIT,TOLD,T,URFT,JS,JN,INTIME,DT,Y,X,XU,YV,GREAT);
     
        % Chama a fun��o PROPS no caso de as propriedades variarem
        if INPRO==1
            [GAMH]=PROPS(NI,NJ,TCON,GAMH);
        end
    
        % Actualiza��o de condi��es de fronteira e fontes se necess�rio
       
        % C�lculo do res�duo normalizado
        SOURCE=(RESORT/SNORM);
    
        % Imprime a informa��o das itera��es no ponto monitor
        fprintf(fid,'%5d',NITER);
        fprintf(fid,'%14.1E',SOURCE);
        fprintf(fid,'%14.3E',T(IMON,JMON));
        fprintf(fid,'%11d',TIME);
        fprintf(fid,'%10d',DT);
        fprintf(fid,'%9d',NSTEP);
        fprintf(fid,'\r \n');
    
        
        % Imprime Temperaturas em intervalos especificados por NITPRI
        if mod(NITER,NITPRI)==0
            PRINT(1,1,NI,NJ,X,Y,T,fid)   

            if NSTEP~=MAXSTP || SOURCE>SORMAX
               % Imprime o r�tulo das informa��es das itera��es no ponto monitor
               fprintf(fid,'NITER        SOURCE     T (%d,%d)    TIME(s)     DT(s)    NSTEP',IMON,JMON);
               fprintf(fid,'\r \n');
            end   
        end
          
        
        % Testa res�duo do processo iterativo
        if SOURCE<SORMAX
            break
        end
        
        % Termina c�culos se a solu��o n�o converge (MAXIT e RES�DUO<10)
        if NITER>=MAXIT && SOURCE>=10
            error ('myApp:argChk','N�o Convergiu segundo o crit�rio especificado \n')
        end
        

    % Termina ciclo no espa�o
    end



    if NSTEP == 3 || NSTEP == 6 || NSTEP == 25 || NSTEP == 100 || NSTEP == 150 || NSTEP == 250
               
        Tmesh(c,1,1) = T(Pts_grid(1,1),Pts_grid(1,2));
        Tmesh(c,2,1) = TIME;
        Tmesh(c,1,2) = T(Pts_grid(2,1),Pts_grid(2,2));
        Tmesh(c,2,2) = TIME;
        Tmesh(c,1,3) = T(Pts_grid(3,1),Pts_grid(3,2));
        Tmesh(c,2,3) = TIME;
        Tmesh(c,1,4) = T(Pts_grid(4,1),Pts_grid(4,2));
        Tmesh(c,2,4) = TIME;
        Tmesh(c,1,5) = T(Pts_grid(5,1),Pts_grid(5,2));
        Tmesh(c,2,5) = TIME;
        Tmesh(c,1,6) = T(Pts_grid(6,1),Pts_grid(6,2));
        Tmesh(c,2,6) = TIME;
        Tmesh(c,1,7) = T(Pts_grid(7,1),Pts_grid(7,2));
        Tmesh(c,2,7) = TIME;
        Tmesh(c,1,8) = T(Pts_grid(8,1),Pts_grid(8,2));
        Tmesh(c,2,8) = TIME;
        Tmesh(c,1,9) = T(Pts_grid(9,1),Pts_grid(9,2));
        Tmesh(c,2,9) = TIME;
        c = c+1;
    end

%     Temperatura_monitor(NSTEP,1) = T(IMON,JMON);    
%     Temperatura_monitor(NSTEP,2) = TIME;
%     Temperatura_centrogeo(NSTEP,1) = T(22,22);
%     Temperatura_centrogeo(NSTEP,2) = TIME;
%     Temperatura_centrosurfN(NSTEP,1) = T(22,42);
%     Temperatura_centrosurfN(NSTEP,2) = TIME;
%     Temperatura_centrosurfO(NSTEP,1) = T(2,22);
%     Temperatura_centrosurfO(NSTEP,2) = TIME;

    
    
    fprintf (fid,'\r \n \r \n \r \n');
    
    % Imprime a solu��o convergida no intervalo especificado por NSTPRI
    if mod (NSTEP,NSTPRI)==0 && mod (NITER,NITPRI)~=0
        PRINT (1,1,NI,NJ,X,Y,T,fid)
    end
       
    if NSTEP~=MAXSTP
       % Imprime o r�tulo das informa��es das itera��es no ponto monitor
       fprintf (fid,'NITER        SOURCE     T (%d,%d)    TIME(s)     DT(s)    NSTEP',IMON,JMON);
       fprintf (fid,'\r \n');
        
    end
    
% Termina ciclo no tempo
end

fclose(fid);

% Desenha gr�fico 2D das Temperaturas--------------------------------------

% Troca os eixos
for jj=1:NI
    for ii=1:NJ
        THI(jj,ii)=T(NI+1-jj,ii);
    end
end

for jj=1:NI
    for ii=1:NJ
        THI2(NJ+1-ii,NI+1-jj)=THI(jj,ii);
    end
    
end
%%% correcao H para W na orientacao x %%%

ZX=W/(NI-1);
ZJ=(H/(NJ-1));

% Plot da Temperatura Final
[X,Y]=meshgrid (0:ZX:W,H:-ZJ:0);

THI2=THI2(2:end-1,2:end-1);
X=X(2:end-1,2:end-1);
Y=Y(2:end-1,2:end-1);

figure(1)
% O n�mero 10 remete para as cores utilizadas
contourf (X,Y,THI2,10);
colorbar;
xlabel('\bfx')
ylabel('\bfy')
zlabel('\bfT')

theta_star_monitor = (Temperatura_monitor(:,1) - 20)/(TINIC-20);
L = 3;
alpha = TCON(1,1)/(DENSIT(1,1)*CV(1,1));
V_corpo = H*W*L;
A_corpo = H*W*2 + H*L*2 + W*L*2;
L_c = 0.12;

Fourier = (alpha*Temperatura_monitor(:,2))/L_c^2;

% figure();
% plot(Fourier,theta_star_monitor)
% ylim([-0.1 1])
% title(['Temperatura adimensionalizada (\theta^*) em fun��o do N�mero de Fourier para T(',num2str(IMON),',',num2str(JMON),')'],'FontSize',14);
% xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
% ylabel('$\theta^*$','Interpreter','latex','FontSize',16)
% 
% 
% theta_star_centrogeo = (Temperatura_centrogeo(:,1) -20)/(TINIC-20);
% 
% figure();
% plot(Fourier,theta_star_centrogeo)
% ylim([-0.1 1])
% title(['Temperatura adimensionalizada (\theta^*) em fun��o do N�mero de Fourier para T(22,22)'],'FontSize',14);
% xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
% ylabel('$\theta^*$','Interpreter','latex','FontSize',16)
% 
% 
% theta_star_centrosurfN = (Temperatura_centrosurfN(:,1) -20)/(TINIC-20);
% 
% figure();
% plot(Fourier,theta_star_centrosurfN)
% ylim([-0.1 1])
% title(['Temperatura adimensionalizada (\theta^*) em fun��o do N�mero de Fourier para T(22,42)'],'FontSize',14);
% xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
% ylabel('$\theta^*$','Interpreter','latex','FontSize',16)
% 
% 
% theta_star_centrosurfO = (Temperatura_centrosurfO(:,1) -20)/(TINIC-20);
% 
% figure();
% plot(Fourier,theta_star_centrosurfO)
% ylim([-0.1 1])
% title(['Temperatura adimensionalizada (\theta^*) em fun��o do N�mero de Fourier para T(2,22)'],'FontSize',14);
% xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
% ylabel('$\theta^*$','Interpreter','latex','FontSize',16)


for s=1:9
    figure()
    plot(Tmesh(:,2,s),Tmesh(:,1,s))
    ylim([0 1150])
    title(['Temperatura em fun��o tempo para T(' ,num2str(Pts_grid(s,1)), ',' ,num2str(Pts_grid(s,2)), ')'],'FontSize',14);
    xlabel('Tempo (s)','Interpreter','latex','FontSize',14)
    ylabel('Temperatura ($^\circ$C)','Interpreter','latex','FontSize',16)
end

%% inven�oes do francisco
% H = 0.24; 
% W = 0.24; 
% L = 3.00; 
% h = 250; %W/m^2K
% T_in = 1150; %�C
% T_amb = 20;
% ro = 7930;
% c = 385;
% k = 121;
% v = 39.6*10^-6;
% alpha = k/(ro*c);
% L_x = 0.12;
% L_y = 0.24;
% L_z = 1.50;
% Bi_x = h*L_x/k;
% Bi_y = h*L_y/k;
% Bi_z = h*L_z/k;
% 
% Bi = [Bi_x Bi_y Bi_z]
% 
% 
% for l=1:3
%     fun = @(csi)csi*tan(csi)-Bi(l);
%     j=1;
%     for i=1:200
%         out = fzero(fun, i-1);
%         if abs(fun(out))<0.05
%             if j==1
%                 csi(l,j)=out;
%                 j = j+1;
%             else if out ~= csi(l,j-1) && out-csi(l,j-1)>0.5
%                     csi(l,j) = out;
%                     j= j+1;
%                 end
%             end
%         end
%     end
% end
% 
% for i=1:20
%     ksi_x(i)=csi(1,i+1);
%     ksi_y(i)=csi(2,i+1);
%     ksi_z(i)=csi(3,i+1);
% end
% x = linspace(0, 1, 6);
% y = linspace(0, 1, 20);
% z = linspace(0, 1, 5);
% theta_estrela_x = 0;
% theta_estrela_y = 0;
% theta_estrela_z = 0;
% 
% V_corpo = H*W*L;
% A_corpo = H*W*2 + H*L*2 + W*L;
% x = [0 1];
% y = [0 0.5 1];
% theta_estrela_x = 0;
% theta_estrela_y = 0;
% theta_estrela_z = 0;
% theta_estrela = 0;
% z = [0 1];
% t = linspace(0, 4000, 200);
% for k=1:length(z)
%     for i=1:2
%         for j=1:3
%             for p=1:length(ksi_x)
%                 C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
%                 C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
%                 C_z(p) = 4*sin(ksi_z(p))/(2*ksi_z(p) + sin(2*ksi_z(p)));
%                 theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
%                 theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(j));
%                 theta_estrela_z = theta_estrela_z + C_z(p)*exp(-ksi_z(p)^2*alpha.*t/(L/2)^2)*cos(ksi_z(p)*z(k));
%             end
%             theta_estrela_2D = theta_estrela_x.*theta_estrela_y;
%             theta_estrela_3D = theta_estrela_x.*theta_estrela_y.*theta_estrela_z;
%             theta_star_lcm = exp(-h*A_corpo/(ro*V_corpo*c).*t);
%             erro_2D_3_D = abs(theta_estrela_3D - theta_estrela_2D);
%             erro_LCM_3D = abs(theta_star_lcm - theta_estrela_3D);
%             theta_estrela_x = 0;
%             theta_estrela_y = 0;
%             theta_estrela_z = 0;
%             figure()
%             if i==1 && j==2
%                 plot(t, theta_estrela_2D,'--',t, theta_estrela_3D, '+', t, theta_star_lcm,'+', Temperatura_monitor(:,2),theta_star_centrogeo)
%                 legend(sprintf('x = %g, y = %g 2D', x(i), y(j)),sprintf('x = %g, y = %g, z = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical')
%                 ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
%                 xlabel("t", 'FontSize', 12)
%             else if i==1 && j==3
%                     plot(t, theta_estrela_2D,'--',t, theta_estrela_3D, '+', t, theta_star_lcm,'+', Temperatura_monitor(:,2),theta_star_centrosurfN)
%                     legend(sprintf('x = %g, y = %g 2D', x(i), y(j)),sprintf('x = %g, y = %g, z = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical')
%                     ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
%                     xlabel("t", 'FontSize', 12)
%                 else if i==2 && j==2
%                         plot(t, theta_estrela_2D,'--',t, theta_estrela_3D, '+', t, theta_star_lcm,'+', Temperatura_monitor(:,2),theta_star_centrosurfO)
%                         legend(sprintf('x = %g, y = %g 2D', x(i), y(j)),sprintf('x = %g, y = %g, z = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical')
%                         ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
%                         xlabel("t", 'FontSize', 12)
%                     else
%                         plot(t, theta_estrela_2D,'--',t, theta_estrela_3D, '+', t, theta_star_lcm,'+')
%                         legend(sprintf('x = %g, y = %g 2D', x(i), y(j)),sprintf('x = %g, y = %g, z = %g 3D', x(i), y(j), z(k)), sprintf('LCM'),'Location','northeast','Orientation','vertical')
%                         ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 18)
%                         xlabel("t", 'FontSize', 12)
%                     end
%                 end
%             end
%         end
%     end
% end


