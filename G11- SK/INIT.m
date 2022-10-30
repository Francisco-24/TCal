%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: INIT (versão MATLAB)
%  Inicializa as matrizes a zero e efectua o cálculo de todos os valores 
%                     associados à malha utilizada
%                                FS
%                             09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[RX,DXEP,DXPW,SEW,XU,RU,RY,DYPS,DYNP,SNS,YV,RV,AN,AS,AE,AW,SU,SP,GAMH,TOLD,T,X,Y]=INIT(INCYLX,INCYLY,NI,NJ,NIM1,NJM1,X,Y,TINIC)

%---------------------- Capitulo 0 - Preliminares ------------------------%


%-------------------- Capitulo 1 - Erro de dados:-------------------------%
%--------------- Programa encerra se INCYLX=1 e INCYLY=1 -----------------%

if INCYLX==1 && INCYLY==1;
    error('myApp:argChk','ERRO-Programa vai encerrar \nINCYLX=INCYLY=1')
    
else
 
    
%--------------- Capitulo 2 - Cálculo dimensões da malha -----------------%

    %Impõe RX=X se axissimétrico na direcção xx   
    for I=1:NI
        RX(I)=1.0;
        if INCYLX==1
           RX(I)=X(I);
        end
    end
    
    %Impõe RY=Y se axissimétrico na direcção yy
    for J=1:NJ
        RY(J)=1.0;
        if INCYLY==1
           RY(J)=Y(J);
        end
    end
    
    %Cálculo da distância entre nos (ver figura 1)
    DXPW(1)=0;
    DXEP(NI)=0;
    for I=1:NIM1
        DXEP(I)=X(I+1)-X(I);
        DXPW(I+1)=DXEP(I);
    end
    
    DYPS(1)=0;
    DYNP(NJ)=0;
    for J=1:NJM1
        DYNP(J)=Y(J+1)-Y(J);
        DYPS(J+1)=DYNP(J);
    end
    
    %Cálculo das dimensões do volume de controlo
    SEW(1)=0;      
    SEW(NI)=0;
    for I=2:NIM1
        SEW(I)=0.5*(DXEP(I)+DXPW(I));  
    end
    
    SNS(1)=0;
    SNS(NJ)=0;
    for J=2:NJM1
        SNS(J)=0.5*(DYNP(J)+DYPS(J));
    end
    
    %Localização das fronteiras dos volumes de controlo
    XU(1)=0;
    RU(1)=0;
    for I=2:NI
        RU(I)=0.5*(RX(I)+RX(I-1));
        XU(I)=0.5*(X(I)+X(I-1));
    end
    
    YV(1)=0;
    RV(1)=0;
    for J=2:NJ
        RV(J)=0.5*(RY(J)+RY(J-1));
        YV(J)=0.5*(Y(J)+Y(J-1));
    end
    
    %Modificação dos valores de fronteira de x e de y
    X(1)=XU(2);
    
    if X(1)<((XU(NI)-XU(2))*1.0E-03)
       X(1)=0;
    end
    
    X(NI)=XU(NI);
    Y(1)=YV(2);
    
    if Y(1)<((YV(NJ)-YV(2))*1.0E-03)
       Y(1)=0.0;
    end   
    
    Y(NJ)=YV(NJ);
    
   
%------------ Capitulo 3 - Inicialização das matrizes a zero -------------%

    J=NJ;
    I=NI;

    AE=zeros(J);
    AW=zeros(J);
    AN=zeros(J);
    AS=zeros(J);
    SU=zeros(J);
    SP=zeros(J);
    GAMH=zeros(I,J);
    TOLD=zeros(I,J);
    
    for ii=1:I
        for jj=1:J
            T(ii,jj)=TINIC;
        end
    
    end
end

%Termina INIT
end 