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
clc

% Utiliza o formato longo (15 dígitos) para maior precisão
format long

% Cria ficheiro para escrever os resultados
fid=fopen ('RESULTS.txt','w');

%---------------------- Capitulo 0 - Preliminares ------------------------%

% Define o número de nós da malha segundo x e y
%%%%% Alterável %%%%%
IT=22;
JT=22;

% Define o número de pontos da malha segundo x e y
%%%%% Alterável %%%%%
NI=22;
NJ=22;

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

Z = zeros(22);
[X,Y] = meshgrid(X,Y)
surface(X,Y,Z)
colormap white


% Estabelece limites do Domínio 
for I=1:NI
    JS (I)=2;
    JN (I)=NJ-1;   
end

% Estabelece ponto monitor segundo x e y (IMON e JMON)
%%%%% Alterável %%%%%
IMON=6;
JMON=6;

% Propriedades do material (Tcond=cond,CV=calor esp,DENSIT=dens), 
% Meio homogéneo
% Aço %%%%% Alterável %%%%%
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

% Parâmetros de controlo do programa
%%%%% Alterável %%%%%
% Número máximo de iterações
MAXIT=10;
% Número máximo de interações no tempo 
MAXSTP=20;
% O output deverá conter os valores de T em intervalos de
NITPRI=110;
% "NITPRI" para "NSTPRI" iterações no tempo
NSTPRI=1;    

% Factor de sub relaxação, Máximo resíduo e intervalo de tempo (s)
%%%%% Alterável %%%%%
URFT=1;
SORMAX=0.001;
DT=50;

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
