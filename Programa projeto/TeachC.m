%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                   Programa TEACH C (versão MATLAB)
%                               CONTROL
%
%     Este programa permite resolver problemas de condução de calor 
%               bidimensional com as seguintes variantes:
%
%    - Coordenadas cartesianas ou coordenadas cilíndricas (axissimétrico)
%    - Regime estacionário ou transiente (não estacionário)
%    - Condutividade térmica (k) uniforme ou variável 
%      (com a temperatura ou com o espaço)
%    - Com a utilização dos três tipos de condições de fronteira existentes
%
%          - Temperatura imposta
%          - Fluxo imposto
%          - Convecção na fronteira do corpo
%
%                                 FS
%                              09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% Apaga variáveis e limpa a consola
clear all
close all
clc

% Utiliza o formato longo (15 dígitos) para maior precisão
format long

% Cria ficheiro para escrever os resultados
fid=fopen ('RESULTS.txt','w');

%---------------------- Capitulo 0 - Preliminares ------------------------%

% Define o número de nós da malha segundo x e y
%%%%% Alterável %%%%%
IT=63;
JT=63;

% Define o número de pontos da malha segundo x e y
%%%%% Alterável %%%%%
NI=63;
NJ=63;

% Constante
GREAT=1.0E30;


%------------ Capitulo 1 - Parâmetros e Índices de controlo --------------%

% Define o tipo de coordenadas (cartesianas- 0 0; 
% cilíndricas vertical 1 0; cilíndricas horizontais 0 1)
%%%%% Alterável %%%%%
INCYLX=0;
INCYLY=0;

% Variável intermédia
NIM1=NI-1;
NJM1=NJ-1;

% Dimensões totais do domínio de solução, largura e altura [m]
%%%%% Alterável %%%%%
W=0.24;
H=0.24;


% Cálculo as abcissas na direção xx
DX(1)=W/(NIM1-1);
X(1)=-0.5*DX(1);
for I=2:NI
    DX(I)=DX(1);
    X(I)=X(I-1)+DX(I-1);
end

% Calcula as abcissas na direção yy
DY(1)=H/(NJM1-1);
Y(1)=-0.5*DY(1);
for J=2:NJ
    DY(J)=DY(1);
    Y(J)=Y(J-1)+DY(J-1);
end

% Estabelece limites do Domínio 
for I=1:NI
    JS (I)=2;
    JN (I)=NJ-1;   
end

% Estabelece ponto monitor segundo x e y (IMON e JMON)
%%%%% Alterável %%%%%
IMON=2;
JMON=62;

% Propriedades do material (Tcond=cond,CV=calor esp,DENSIT=dens), 
% Meio homogéneo
% Aço %%%%% Alterável %%%%%
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

% Parâmetros de controlo do programa
%%%%% Alterável %%%%%
% Número máximo de iterações
MAXIT=400;
% Número máximo de interações no tempo 
MAXSTP=300;
% O output deverá conter os valores de T em intervalos de
NITPRI=400;
% "NITPRI" para "NSTPRI" iterações no tempo
NSTPRI=1;    

% Factor de sub relaxação, Máximo resíduo e intervalo de tempo (s)
%%%%% Alterável %%%%%
URFT=1;
SORMAX=0.001;
DT=20;

% Selecciona o Regime---Estacionário->INTIME=0, Transiente->INTIME=1
%%%%% Alterável %%%%%
INTIME=1;

if INTIME==0
    MAXSTP=1;
end

% Indica se as propriedades são constantes --- constantes->INPRO=0,
% variáveis->INPRO=1
%%%%% Alterável %%%%%
INPRO=0;


%------------------ Capitulo 2 - Opera\[CCedilla]\[OTilde]es iniciais ----------------------%

% Calcula dimensões da malha e anula vectores/matrizes
% Chama função INIT
[RX,DXEP,DXPW,SEW,XU,RU,RY,DYPS,DYNP,SNS,YV,RV,AN,AS,AE,AW,SU,SP,GAMH,TOLD,T,X,Y]=INIT (INCYLX,INCYLY,NI,NJ,NIM1,NJM1,X,Y,TINIC);
TIME=0.0;
 
% Impõe valores de fronteira e inicializa variável dependente

% Valores de fronteira
%%%%% Alterável %%%%%
TTOP=20;
TBOT=1150;
TLEFT=0;
TRIGHT=0;

% Uma vez que de acordo com o enunciado não há temperatura imposta em
% nenhuma das faces estes valores não irão ter importância nos resultados,
% estando apenas a adicionar tempo de computação

% for I=2:NIM1
%     T (I,1)=TBOT;
%     T (I,NJ)=TTOP;
% end
%  
% for J=2:NJM1
%     T (1,J)=TLEFT;
%     T (NI,J)=TRIGHT;
% end
 

% Inicializa variável dependente
% Inicializa campo de propriedades do material
% Chama função PROPS
[GAMH]=PROPS(NI,NJ,TCON,GAMH);
 
 
% Cálculo do factor de normalização do resíduo
AK=BK;
SNORM=AK*(TTOP-TBOT)*W/H;
SNORM=abs(SNORM);

% Escreve as especificações do problema
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

% Chama a função PRINT para imprimir o campo de temperaturas inicial
PRINT(1,1,NI,NJ,X,Y,T,fid)


%-------------- Capitulo 3 - Iterações no tempo e no espaço --------------%

% Indica os pontos de controlo
% Imprime o rótulo das informações das iterações no ponto monitor
fprintf (fid,'NITER        SOURCE        T (%d,%d)    TIME(s)     DT(s)    NSTEP',IMON,JMON);
fprintf (fid,'\r \n');

Temperatura_monitor = zeros(MAXSTP,2);
Temperatura_centrogeo = zeros(MAXSTP,2);
Temperatura_centrosurfN = zeros(MAXSTP,2);
Temperatura_centrosurfO = zeros(MAXSTP,2);
Temperatura_centrosurfS = zeros(MAXSTP,2);
%Tmesh = zeros(8,2,9);
%c = 1;

%Para malha de 82x82
% Pts_grid = [5,78;6,77;5,38;6,37;5,6;6,5;21,78;22,77;21,38;22,37;21,6;22,5;37,78;38,77;37,38;38,37;37,6;38,5];

%Para malha de 62x62
% Pts_grid = [4,59;5,58;4,29;5,28;4,5;5,4;16,59;17,58;16,29;17,28;16,5;17,4;28,59;29,58;28,29;29,28;28,5;29,4];

%Para malha de 42x42
% Pts_grid = [3,40;4,39;3,20;4,19;3,4;4,3;11,40;12,39;11,20;12,19;11,4;12,3;19,40;20,39;19,20;20,19;19,4;20,3];

%Para malha de 22x22
% Pts_grid = [2,21;3,20;2,11;3,10;2,3;3,2;6,21;7,20;6,11;7,10;6,2;7,3;10,21;11,20;10,11;11,10;10,2;11,3];

%Para malha de 12x12
%Pts_grid = [2,11;2,6;2,2;4,11;4,6;4,2;6,11;6,6;6,2];

% Iterações no tempo
for NSTEP=1:MAXSTP
    TIME=TIME+DT;
    for I=1:NI
        for J=1:NJ
            TOLD(I,J)=T(I,J);
        end
    end
  
    % Iterações no espaço
    for NITER=1:MAXIT
        
        % Chama a função CALCT para o cálculo das temperaturas
        [AN,AS,AE,AW,SU,SP,GAMH,CV,DENSIT,TOLD,T,RESORT]=CALCT(NI,NJ,NIM1,NJM1,RX,DXEP,SEW,RU,RY,DYNP,SNS,RV,AN,AS,AE,AW,SU,SP,GAMH,CV,DENSIT,TOLD,T,URFT,JS,JN,INTIME,DT,Y,X,XU,YV,GREAT);
     
        % Chama a função PROPS no caso de as propriedades variarem
        if INPRO==1
            [GAMH]=PROPS(NI,NJ,TCON,GAMH);
        end
    
        % Actualização de condições de fronteira e fontes se necessário
       
        % Cálculo do resíduo normalizado
        SOURCE=(RESORT/SNORM);
    
        % Imprime a informação das iterações no ponto monitor
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
               % Imprime o rótulo das informações das iterações no ponto monitor
               fprintf(fid,'NITER        SOURCE     T (%d,%d)    TIME(s)     DT(s)    NSTEP',IMON,JMON);
               fprintf(fid,'\r \n');
            end   
        end
          
        
        % Testa resíduo do processo iterativo
        if SOURCE<SORMAX
            break
        end
        
        % Termina cáculos se a solução não converge (MAXIT e RESÍDUO<10)
        if NITER>=MAXIT && SOURCE>=10
            error ('myApp:argChk','Não Convergiu segundo o critério especificado \n')
        end
        

    % Termina ciclo no espaço
    end


%Para malha 82x82, 62x62, 42x42 e 22x22
%     if NSTEP == 3 || NSTEP == 6 || NSTEP == 25 || NSTEP == 50 || NSTEP == 100 || NSTEP == 150 || NSTEP == 200 || NSTEP == 250
%                
%         Tmesh(c,1,1) = (T(Pts_grid(1,1),Pts_grid(1,2))+T(Pts_grid(2,1),Pts_grid(2,2)))/2;
%         Tmesh(c,2,1) = TIME;
%         Tmesh(c,1,2) = (T(Pts_grid(3,1),Pts_grid(3,2))+T(Pts_grid(4,1),Pts_grid(4,2)))/2;
%         Tmesh(c,2,2) = TIME;
%         Tmesh(c,1,3) = (T(Pts_grid(5,1),Pts_grid(5,2))+T(Pts_grid(6,1),Pts_grid(6,2)))/2;
%         Tmesh(c,2,3) = TIME;
%         Tmesh(c,1,4) = (T(Pts_grid(7,1),Pts_grid(7,2))+T(Pts_grid(8,1),Pts_grid(8,2)))/2;
%         Tmesh(c,2,4) = TIME;
%         Tmesh(c,1,5) = (T(Pts_grid(9,1),Pts_grid(9,2))+T(Pts_grid(10,1),Pts_grid(10,2)))/2;
%         Tmesh(c,2,5) = TIME;
%         Tmesh(c,1,6) = (T(Pts_grid(11,1),Pts_grid(11,2))+T(Pts_grid(12,1),Pts_grid(12,2)))/2;
%         Tmesh(c,2,6) = TIME;
%         Tmesh(c,1,7) = (T(Pts_grid(13,1),Pts_grid(13,2))+T(Pts_grid(14,1),Pts_grid(14,2)))/2;
%         Tmesh(c,2,7) = TIME;
%         Tmesh(c,1,8) = (T(Pts_grid(15,1),Pts_grid(15,2))+T(Pts_grid(16,1),Pts_grid(16,2)))/2;
%         Tmesh(c,2,8) = TIME;
%         Tmesh(c,1,9) = (T(Pts_grid(17,1),Pts_grid(17,2))+T(Pts_grid(18,1),Pts_grid(18,2)))/2;
%         Tmesh(c,2,9) = TIME;
%         c = c+1;
%     end
    
%Para malha 12x12
%     if NSTEP == 3 || NSTEP == 6 || NSTEP == 25 || NSTEP == 100 || NSTEP == 150 || NSTEP == 250
%                
%         Tmesh(c,1,1) = T(Pts_grid(1,1),Pts_grid(1,2));
%         Tmesh(c,2,1) = TIME;
%         Tmesh(c,1,2) = T(Pts_grid(2,1),Pts_grid(2,2));
%         Tmesh(c,2,2) = TIME;
%         Tmesh(c,1,3) = T(Pts_grid(3,1),Pts_grid(3,2));
%         Tmesh(c,2,3) = TIME;
%         Tmesh(c,1,4) = T(Pts_grid(4,1),Pts_grid(4,2));
%         Tmesh(c,2,4) = TIME;
%         Tmesh(c,1,5) = T(Pts_grid(5,1),Pts_grid(5,2));
%         Tmesh(c,2,5) = TIME;
%         Tmesh(c,1,6) = T(Pts_grid(6,1),Pts_grid(6,2));
%         Tmesh(c,2,6) = TIME;
%         Tmesh(c,1,7) = T(Pts_grid(7,1),Pts_grid(7,2));
%         Tmesh(c,2,7) = TIME;
%         Tmesh(c,1,8) = T(Pts_grid(8,1),Pts_grid(8,2));
%         Tmesh(c,2,8) = TIME;
%         Tmesh(c,1,9) = T(Pts_grid(9,1),Pts_grid(9,2));
%         Tmesh(c,2,9) = TIME;
%         c = c+1;
%     end

    Temperatura_monitor(NSTEP,1) = T(IMON,JMON);    
    Temperatura_monitor(NSTEP,2) = TIME;
    Temperatura_centrogeo(NSTEP,1) = T(32,32);
    Temperatura_centrogeo(NSTEP,2) = TIME;
    Temperatura_centrosurfN(NSTEP,1) = T(32,62);
    Temperatura_centrosurfN(NSTEP,2) = TIME;
    Temperatura_centrosurfO(NSTEP,1) = T(2,32);
    Temperatura_centrosurfO(NSTEP,2) = TIME;
    Temperatura_centrosurfS(NSTEP,1) = T(32,2);
    Temperatura_centrosurfS(NSTEP,2) = TIME;

    
    
    fprintf (fid,'\r \n \r \n \r \n');
    
    % Imprime a solução convergida no intervalo especificado por NSTPRI
    if mod (NSTEP,NSTPRI)==0 && mod (NITER,NITPRI)~=0
        PRINT (1,1,NI,NJ,X,Y,T,fid)
    end
       
    if NSTEP~=MAXSTP
       % Imprime o rótulo das informações das iteraÇÕes no ponto monitor
       fprintf (fid,'NITER        SOURCE     T (%d,%d)    TIME(s)     DT(s)    NSTEP',IMON,JMON);
       fprintf (fid,'\r \n');
        
    end
    
% Termina ciclo no tempo
end

% for sheet = 1:9
%   filename = 'Grid_independency.xlsx';
%   
%   Para malha 82x82
%   writematrix(Tmesh(:,:,sheet), filename,'Sheet',sheet,'Range', 'P4:Q11')
% 
%   Para malha 62x62
%   writematrix(Tmesh(:,:,sheet), filename,'Sheet',sheet,'Range', 'M4:N11')
% 
%   Para malha 42x42
%   writematrix(Tmesh(:,:,sheet), filename,'Sheet',sheet,'Range', 'J4:K11')
% 
%   Para malha 22x22
%   writematrix(Tmesh(:,:,sheet), filename,'Sheet',sheet,'Range', 'G4:H11')
% 
%   Para malha 12x12
%   writematrix(Tmesh(:,:,sheet), filename,'Sheet',sheet,'Range', 'D4:E11')
% 
% end

fclose(fid);

% Desenha gráfico 2D das Temperaturas--------------------------------------

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
% O número 10 remete para as cores utilizadas
contourf (X,Y,THI2,10);
colorbar;
xlabel('\bfx')
ylabel('\bfy')
zlabel('\bfT')



%Cálculo do número de Fourier
L = 3;
alpha = TCON(1,1)/(DENSIT(1,1)*CV(1,1));
L_c = 0.12;
Fourier = (alpha*Temperatura_monitor(:,2))/L_c^2;

%Plots

theta_star_centrogeo = (Temperatura_centrogeo(:,1) -20)/(TINIC-20);

figure();
plot(Fourier,theta_star_centrogeo)
title(['Temperatura adimensionalizada (\theta^*) em função do Número de Fourier para T(32,32)'],'FontSize',14);
xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
ylabel('$\theta^*$','Interpreter','latex','FontSize',16)


theta_star_centrosurfN = (Temperatura_centrosurfN(:,1) -20)/(TINIC-20);

figure();
plot(Fourier,theta_star_centrosurfN)
title(['Temperatura adimensionalizada (\theta^*) em função do Número de Fourier para T(32,62)'],'FontSize',14);
xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
ylabel('$\theta^*$','Interpreter','latex','FontSize',16)


theta_star_centrosurfO = (Temperatura_centrosurfO(:,1) -20)/(TINIC-20);

figure();
plot(Fourier,theta_star_centrosurfO)
title(['Temperatura adimensionalizada (\theta^*) em função do Número de Fourier para T(2,32)'],'FontSize',14);
xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
ylabel('$\theta^*$','Interpreter','latex','FontSize',16)

theta_star_centrosurfS = (Temperatura_centrosurfS(:,1) -20)/(TINIC-20);

figure();
plot(Fourier,theta_star_centrosurfS)
title(['Temperatura adimensionalizada (\theta^*) em função do Número de Fourier para T(32,2)'],'FontSize',14);
xlabel('N{\''{u}}mero de Fourier','Interpreter','latex','FontSize',14)
ylabel('$\theta^*$','Interpreter','latex','FontSize',16)

%% Método analitico && comparação
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


for i=1:20
    ksi_x(i)=csi(1,i+1);
    ksi_y(i)=csi(2,i+1);
    ksi_z(i)=csi(3,i+1);
end
x = linspace(0, 1, 6);
y = linspace(0, 1, 20);
z = linspace(0, 1, 5);
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;

V_corpo = H*W*L;
A_corpo = H*W*2 + H*L*2 + W*L;
x = [0 1];
y = [0 0.5 1];
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;
theta_estrela = 0;
z = [0 1];
t = linspace(0, 6000, 301);
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

            if i==1 && j==1 && k==1
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrosurfS,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==1 && j==1 && k==2
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-',t*alpha/(H/2)^2, theta_estrela_3D, 'm--', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrosurfS,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)),sprintf('x* = %g, y* = %g, z* = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==1 && j==2 && k==1
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrogeo,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==1 && j==2 && k==2
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-',t*alpha/(H/2)^2, theta_estrela_3D, 'm--', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrogeo,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)),sprintf('x* = %g, y* = %g, z* = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==1 && j==3 && k==1
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrosurfN,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==1 && j==3 && k==2
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-',t*alpha/(H/2)^2, theta_estrela_3D, 'm--', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrosurfN,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)),sprintf('x* = %g, y* = %g, z* = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==2 && j==2 && k==1
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrosurfO,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            elseif i==2 && j==2 && k==2
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-',t*alpha/(H/2)^2, theta_estrela_3D, 'm--', t*alpha/(H/2)^2, theta_star_lcm,'y-', alpha*Temperatura_monitor(:,2)/L_c^2,theta_star_centrosurfO,'b--','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)),sprintf('x* = %g, y* = %g, z* = %g 3D', x(i), y(j), z(k)), sprintf('LCM'), 'Teach C','Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            else
                plot(t*alpha/(H/2)^2, theta_estrela_2D,'r-',t*alpha/(H/2)^2, theta_estrela_3D, 'm--', t*alpha/(H/2)^2, theta_star_lcm,'y-','LineWidth',1.5)
                legend(sprintf('x* = %g, y* = %g 2D', x(i), y(j)),sprintf('x* = %g, y* = %g, z* = %g 3D', x(i), y(j), z(k)), sprintf('LCM'),'Location','northeast','Orientation','vertical','FontSize', 15)
                ylabel("$\theta*$", 'Interpreter','latex', 'FontSize', 20)
                xlabel("Fo", 'FontSize', 20)
            end
        end
    end
end


% cálculo dos erros 
x = linspace(0, 1, 6);
y = linspace(0, 1, 20);
z = linspace(0, 1, 5);
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;

V_corpo = H*W*L;
A_corpo = H*W*2 + H*L*2 + W*L;
x = [0 1];
y = [0 0.5 1];
theta_estrela_x = 0;
theta_estrela_y = 0;
theta_estrela_z = 0;
theta_estrela = 0;
z = 0;
t = linspace(20, 6000, 300);

for i=1:2
    for j=1:3
        for p=1:length(ksi_x)
            C_x(p) = 4*sin(ksi_x(p))/(2*ksi_x(p) + sin(2*ksi_x(p)));
            C_y(p) = 4*sin(ksi_y(p))/(2*ksi_y(p) + sin(2*ksi_y(p)));
            C_z(p) = 4*sin(ksi_z(p))/(2*ksi_z(p) + sin(2*ksi_z(p)));
            theta_estrela_x = theta_estrela_x + C_x(p)*exp(-ksi_x(p)^2*alpha.*t/(H/2)^2)*cos(ksi_x(p)*x(i));
            theta_estrela_y = theta_estrela_y + C_y(p)*exp(-ksi_y(p)^2*alpha.*t/(W)^2)*cos(ksi_y(p)*y(j));
            theta_estrela_z = theta_estrela_z + C_z(p)*exp(-ksi_z(p)^2*alpha.*t/(L/2)^2)*cos(ksi_z(p)*z);
        end

        theta_estrela_2D = theta_estrela_x.*theta_estrela_y;
        theta_estrela_3D = theta_estrela_x.*theta_estrela_y.*theta_estrela_z;
        theta_estrela_2D_t = theta_estrela_2D.';
        theta_estrela_3D = theta_estrela_3D.';
        theta_star_lcm = exp(-h*A_corpo/(ro*V_corpo*c).*t);
        erro_2D_3_D = abs(theta_estrela_3D - theta_estrela_2D);
        erro_LCM_3D = abs(theta_star_lcm - theta_estrela_2D);
        theta_estrela_x = 0;
        theta_estrela_y = 0;
        theta_estrela_z = 0;
        
        if i==1 && j==1
            erro_teachC = abs(theta_star_centrosurfS-theta_estrela_2D_t);
            figure()
            plot(t*alpha/(H/2)^2, erro_teachC,'LineWidth',1.5)
            title('Erro absoluto (solução numérica vs solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
            
            figure()
            plot(t*alpha/(H/2)^2, erro_LCM_3D,'LineWidth',1.5)
            title('Erro absoluto (método da Capacitância Global vs Solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
        elseif i==1 && j==2
            erro_teachC = abs(theta_star_centrogeo-theta_estrela_2D_t);
            figure()
            plot(t*alpha/(H/2)^2, erro_teachC,'LineWidth',1.5)
            title('Erro absoluto (solução numérica vs solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
            
            figure()
            plot(t*alpha/(H/2)^2, erro_LCM_3D,'LineWidth',1.5)
            title('Erro absoluto (método da Capacitância Global vs Solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
        elseif i==1 && j==3
            erro_teachC = abs(theta_star_centrosurfN-theta_estrela_2D_t);
            figure()
            plot(t*alpha/(H/2)^2, erro_teachC,'LineWidth',1.5)
            title('Erro absoluto (solução numérica vs solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
            
            figure()
            plot(t*alpha/(H/2)^2, erro_LCM_3D,'LineWidth',1.5)
            title('Erro absoluto (método da Capacitância Global vs Solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
        elseif i==2 && j==2
            erro_teachC = abs(theta_star_centrosurfO-theta_estrela_2D_t);
            figure()
            plot(t*alpha/(H/2)^2, erro_teachC,'LineWidth',1.5)
            title('Erro absoluto (solução numérica vs solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
            
            figure()
            plot(t*alpha/(H/2)^2, erro_LCM_3D,'LineWidth',1.5)
            title('Erro absoluto (método da Capacitância Global vs Solução Analítica)', 'FontSize', 20)
            legend(sprintf('x* = %g, y* = %g', x(i), y(j)),'Location','northeast','Orientation','vertical','FontSize', 10)
            ylabel("erro absoluto", 'Interpreter','latex', 'FontSize', 20)
            xlabel("Fo", 'FontSize', 20)
        end
    end
end

