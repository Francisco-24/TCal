%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: SOLVE (vers�o MATLAB)
%  Esta fun��o tem por finalidade a resolu��o de um sistema de equa��es em
%        que a matriz dos coeficientes � uma matriz tri-diagonal
%                                FS
%                             09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[PHI]=SOLVE(NI,NJ,IL,JSTART,JEND,PHI,AN,AS,AE,AW,AP,SU,SP)

%////////////////////// Algoritmo de THOMA: TDMA \\\\\\\\\\\\\\\\\\\\\\\\\%


%--------------------- Capitulo 0 - Preliminares -------------------------%


%------------ Capitulo 1 - Processo de itera�ao linha a linha ------------%

JENDM1=JEND-1;
JSTM1=JSTART-1;

%Inicializa��o dos coeficientes
for J=JSTART:JEND
    A(J)=0;
    B(J)=0;
    C(J)=0;
    D(J)=0;
    
end

A(JSTM1)=0;
I=IL;
C(JSTM1)=PHI(I,JSTM1);

for J=JSTART:JENDM1
    
    %Assemblagem dos coeficientes TDMA
    A(J)=AN(J);
    B(J)=AS(J);
    C(J)=AE(J)*PHI(I+1,J)+AW(J)*PHI(I-1,J)+SU(J);
    D(J)=AP(J);
    
    %Calcula os coeficientes da f�rmula de recorr�ncia
    DIV=D(J)-B(J)*A(J-1);
    TERM=1.0/DIV;
    A(J)=A(J)*TERM;
    C(J)=(C(J)+B(J)*C(J-1))*TERM;
    
end

%Obtem novos PHI's por substitui��o para tr�s
for JJ=JSTART:JENDM1
    J=JEND+JSTM1-JJ;
    PHI(I,J)=A(J)*PHI(I,J+1)+C(J);
    
end


%Termina SOLVE
end