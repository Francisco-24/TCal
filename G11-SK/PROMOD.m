%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: PROMOD (versão MATLAB)
%                Define e introduz as Condições de Fronteira 
%                    -Temperatura Imposta
%                    -Fluxo Imposto
%                    -Convecção
%                    -Simetria
%                                FS
%                             09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[AN,AS,AE,AW,SU,SP,T]=PROMOD(NI,NJ,NIM1,NJM1,IL,RV,YV,Y,SNS,SEW,X,XU,AN,AS,AE,AW,SU,SP,GAMH,T,RY)


%----------------------- Capitulo 0 - Preliminares -----------------------%


%------------------ Capitulo 1 - Condição de Fronteira -------------------%

%O Utilizador define o tipo de condição de fronteira (1-Temperatura, 
%2-Fluxo, 3-Convecção e 4-simetria para a fronteira Norte, Sul, Este e 
%Oeste (CFN,CFS,CFW,CFE)
%%%%% Alterável %%%%%
CFN=3;
CFS=4;
CFW=3;
CFE=3;

%Identifica o tipo de condição de fronteira
I=IL;

%%%FRONTEIRA NORTE---------------------------------------------------------
%Temperatura imposta
if CFN==1
       
    RDYN=RV(NJ)/(YV(NJ)-Y(NJM1));
    AN(NJM1)=0;
    DN=GAMH(IL,NJM1)*SEW(IL)*RDYN;
    SU(NJM1)=SU(NJM1)+DN*T(IL,NJ);
    SP(NJM1)=SP(NJM1)-DN;

    
%Fluxo imposto
elseif CFN==2
    
    QN=20;           %%%%% alterável %%%%%
    AN(NJM1)=0;
    DN=QN*SEW(IL)*RV(NJ);
    SU(NJM1)=SU(NJM1)+DN;
    SP(NJM1)=SP(NJM1);

    
%Convecção
elseif CFN==3
    
    HCONV=250;     %%%%% alterável %%%%%
    TF=20;          %%%%% alterável %%%%%
    RDYN=YV(NJ)-Y(NJM1);
    AN(NJM1)=0;
    DN1=RDYN/GAMH(IL,NJM1)+1/HCONV;
    DN=SEW(IL)*RV(NJ)/DN1;
    SU(NJM1)=SU(NJM1)+DN*TF;
    SP(NJM1)=SP(NJM1)-DN;

    
%simetria
elseif CFN==4
    
    AN(NJM1)=0;
    T(IL,NJ)=T(IL,NJM1);
 
    
%Caso o utilizador tenha introduzido incorrectamente a condição de
%fronteira o programa termina
else
    error('myApp:argChk','Condição de fronteira mal introduzida, o programa vai encerrar')

    
end
%%%------------------------------------------------------------------------

%%%FRONTEIRA SUL-----------------------------------------------------------
%Temperatura imposta
if CFS==1
    
    RDYS=RV(2)/(Y(2)-YV(2));
    AS(2)=0;
    DS=GAMH(IL,2)*SEW(IL)*RDYS;
    SU(2)=SU(2)+DS*T(IL,1);
    SP(2)=SP(2)-DS;

    
%Fluxo imposto
elseif CFS==2

    QS=0;          %%%%% alterável %%%%%
    AS(2)=0;
    DS=QS*SEW(IL)*RV(2);
    SU(2)=SU(2)+DS;
    SP(2)=SP(2);
    
elseif CFS==3
    
    HCONV=12.5;     %%%%% alterável %%%%%
    TF=80;          %%%%% alterável %%%%%
    RDYS=Y(2)-YV(2);
    AS(2)=0;
    DS1=RDYS/GAMH(IL,2)+1/HCONV;
    DS=SEW(IL)*RV(2)/DS1;
    SU(2)=SU(2)+DS*TF;
    SP(2)=SP(2)-DS;
    

%simetria
elseif CFS==4
    
    AS(2)=0;
    T(IL,2)=T(IL,2);
    
    
%Caso o utilizador tenha introduzido incorrectamente a condição de
%fronteira o programa termina    
else
    
    error('myApp:argChk','Condição de fronteira mal introduzida, o programa vai encerrar')
   
    
end
%%%------------------------------------------------------------------------

%%%FRONTEIRA OESTE---------------------------------------------------------
%Temperatura imposta
if CFW==1
    
    if IL==2
        DXW=X(2)-XU(2);
        for J=2:NJM1
            AW(J)=0;
            DW=GAMH(IL,J)*SNS(J)*RY(J)/DXW;
            SU(J)=SU(J)+DW*T(1,J);
            SP(J)=SP(J)-DW;
        end
    end

    
%Fluxo imposto
elseif CFW==2
    
    if IL==2
        QW=20;           %%%%% alterável %%%%%       
        for J=2:NJM1
            AW(J)=0;
            DW=QW*SNS(J)*RY(J);   
            SU(J)=SU(J)+DW;
            SP(J)=SP(J);
            
        end
    end
  
    
%Convecçâo
elseif CFW==3
    
     if IL==2
        HCONV=250;     %%%%% alterável %%%%%
        TF=20;          %%%%% alterável %%%%%
        DXW=X(2)-XU(2);
        for J=2:NJM1
            AW(J)=0;
            DW1=DXW/GAMH(IL,J)+1/HCONV;
            DW=SNS(J)*RY(J)/DW1;
            SU(J)=SU(J)+DW*TF;
            SP(J)=SP(J)-DW;
        end
     end
    

elseif CFW==4
    
     if IL==2
        
        for J=2:NJM1
        AW(J)=0;
        T(1,J)=T(2,J);
        end
        
     end
    
     
%Caso o utilizador tenha introduzido incorrectamente a condição de
%fronteira o programa termina   
else
    
    error('myApp:argChk','Condição de fronteira mal introduzida, o programa vai encerrar')
  
    
end
%%%------------------------------------------------------------------------    
 
%%%FRONTEIRA ESTE----------------------------------------------------------
%Temperatura imposta
if CFE==1
    
    if IL==NIM1
        DXE=XU(NI)-X(NIM1);
        for J=2:NJM1
            AE(J)=0;
            DE=GAMH(IL,J)*SNS(J)*RY(J)/DXE;
            SU(J)=SU(J)+DE*T(NI,J);
            SP(J)=SP(J)-DE;
            
        end  
    end
    
    
%Fluxo imposto
elseif CFE==2
    
     if IL==NIM1
        QE=0;            %%%%% alterável %%%%%
        for J=2:NJM1
            AE(J)=0;
            DE=QE*SNS(J)*RY(J);  
            SU(J)=SU(J)+DE;
            SP(J)=SP(J);
            
        end
     end 
    
     
%Convecção
elseif CFE==3
    
    if IL==NIM1
        HCONV=250;     %%%%% alterável %%%%%
        TF=20;          %%%%% alterável %%%%%
        DXE=XU(NI)-X(NIM1);
        for J=2:NJM1
            AE(J)=0;
            DE1=DXE/GAMH(IL,J)+1/HCONV;
            DE=SNS(J)*RY(J)/DE1;
            SU(J)=SU(J)+DE*TF;
            SP(J)=SP(J)-DE;
        end
    end
    
     
%simetria    
elseif CFE==4
    
     if IL==NIM1
        
        for J=2:NJM1
        AE(J)=0;
        T(NI,J)=T(NIM1,J);
        end
        
     end
    

%Caso o utilizador tenha introduzido incorrectamente a condição de
%fronteira o programa termina   
else
     
    error('myApp:argChk','Condição de fronteira mal introduzida, o programa vai encerrar')
 
    
end
%%%------------------------------------------------------------------------


%Termina PROMOD
end