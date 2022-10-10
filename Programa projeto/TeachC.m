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
clc

% Utiliza o formato longo (15 d�gitos) para maior precis�o
format long

% Cria ficheiro para escrever os resultados
fid=fopen ('RESULTS.txt','w');

%---------------------- Capitulo 0 - Preliminares ------------------------%

% Define o n�mero de n�s da malha segundo x e y
%%%%% Alter�vel %%%%%
IT=22;
JT=22;

% Define o n�mero de pontos da malha segundo x e y
%%%%% Alter�vel %%%%%
NI=22;
NJ=22;

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

Z = zeros(22);
[X,Y] = meshgrid(X,Y)
surface(X,Y,Z)
colormap white


% Estabelece limites do Dom�nio 
for I=1:NI
    JS (I)=2;
    JN (I)=NJ-1;   
end

% Estabelece ponto monitor segundo x e y (IMON e JMON)
%%%%% Alter�vel %%%%%
IMON=6;
JMON=6;

% Propriedades do material (Tcond=cond,CV=calor esp,DENSIT=dens), 
% Meio homog�neo
% A�o %%%%% Alter�vel %%%%%
for I=1:NI
    for J=1:NJ
        TCON (I,J)=14.68;
        CV (I,J)=485.67;
        DENSIT (I,J)=7800;
        if I==1 && J==1
           BK=TCON (I,J);
            
        end
    end
end

% Temperatura inicial
TINIC=0;

% Par�metros de controlo do programa
%%%%% Alter�vel %%%%%
% N�mero m�ximo de itera��es
MAXIT=10;
% N�mero m�ximo de intera��es no tempo 
MAXSTP=20;
% O output dever� conter os valores de T em intervalos de
NITPRI=110;
% "NITPRI" para "NSTPRI" itera��es no tempo
NSTPRI=1;    

% Factor de sub relaxa��o, M�ximo res�duo e intervalo de tempo (s)
%%%%% Alter�vel %%%%%
URFT=1;
SORMAX=0.001;
DT=50;

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
TTOP=100;
TBOT=0;
TLEFT=0;
TRIGHT=0;
 
for I=2:NIM1
    T (I,1)=TBOT;
    T (I,NJ)=TTOP;
end
 
for J=2:NJM1
    T (1,J)=TLEFT;
    T (NI,J)=TRIGHT;
end
 
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
