%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: CALCT (vers�o MATLAB)
%      A finalidade deste bloco � a de calcular os coeficientes da 
%               equa��o diferencial de calor discretizada
%                                FS
%                             09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[AN,AS,AE,AW,SU,SP,GAMH,CV,DENSIT,TOLD,T,RESORT]=CALCT(NI,NJ,NIM1,NJM1,RX,DXEP,SEW,RU,RY,DYNP,SNS,RV,AN,AS,AE,AW,SU,SP,GAMH,CV,DENSIT,TOLD,T,URFT,JS,JN,INTIME,DT,Y,X,XU,YV,GREAT)
                                                         
%--------------------- Capitulo 0 - Preliminares -------------------------%


%------------------ Capitulo 1 - Calculo coeficientes --------------------%

RESORT=0;

%Ciclo pelas colunas I=cte
for I=2:NIM1

    %Encontra limites JJ inferior e superior para cada coluna I=cte
    LJS=JS(I);
    LJN=JN(I);
    
    %Calcula coeficientes para toda a coluna I=cte
    for J=LJS:LJN
        
        %Determina �reas e volume
        AREAN=RV(J+1)*SEW(I)*RX(I);
        AREAE=RY(J)*SNS(J)*RU(I+1);
        VOL=RY(J)*SNS(J)*SEW(I)*RX(I);
        
        %Calcula coeficientes de difus�o
        GAMN=0.5*(GAMH(I,J)+GAMH(I,J+1));
        GAME=0.5*(GAMH(I,J)+GAMH(I+1,J));
        DN=GAMN*AREAN/DYNP(J);
        DE=GAME*AREAE/DXEP(I);
        
        %Termos de fonte quando existentes
        %%%%% Alter�vel %%%%%
        SU(J)=0;
        SP(J)=0;
        
        %Coeficientes Transientes
        if INTIME==1
            DP=VOL*CV(I,J)*DENSIT(I,J)/DT;
            SU(J)=SU(J)+DP*TOLD(I,J);
            SP(J)=SP(J)-DP;
        end
        
        %Calcula coeficientes
        AN(J)=DN;
        AS(J)=AN(J-1);
        AW(J)=AE(J);
        AE(J)=DE;
        
    end
    
    IL=I;
     
    
%----- Capitulo 2 - Modifica��es do problema: condi��es de fronteira -----%

    %Chama a Fun��o PROMOD
    [AN,AS,AE,AW,SU,SP,T]=PROMOD(NI,NJ,NIM1,NJM1,IL,RV,YV,Y,SNS,SEW,X,XU,AN,AS,AE,AW,SU,SP,GAMH,T,RY);
    
    
%-------------- Capitulo 3 - Coeficientes finais + Res�duo ---------------%

    for J=LJS:LJN
        AP(J)=AN(J)+AS(J)+AE(J)+AW(J)-SP(J);
        RESOR=AN(J)*T(I,J+1)+AS(J)*T(I,J-1)+AE(J)*T(I+1,J)+AW(J)*T(I-1,J)-AP(J)*T(I,J)+SU(J);
        VOL=RY(J)*SEW(I)*SNS(J)*RX(I);
        
        %Modifica RESOR se condi��es de fronteira s�o aplicadas com recurso
        %a SP=-GREAT
        if -SP(J)>(0.5*GREAT)
            RESOR=RESOR/GREAT;
        end
        
        %Soma res�duo para esta coluna I=cte
        RESORT=RESORT+abs(RESOR);
        
        %Sub-relaxa��o
        AP(J)=AP(J)/URFT;
        SU(J)=SU(J)+(1.0-URFT)*AP(J)*T(I,J);
    end
    
    
%------------- Capitulo 4 - Solu��o das equa��es alg�bricas --------------%

    %Efectua itera��o linha-a-linha
    [T]=SOLVE(NI,NJ,IL,LJS,LJN+1,T,AN,AS,AE,AW,AP,SU,SP);   
               
end

%Termina CALCT
end