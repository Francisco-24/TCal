%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: PRINT (versão MATLAB)
%                            Produz o Output 
%                     Mostra ao utilizador a matriz T/PHI 
%                                FS
%                             09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[]=PRINT(ISTART,JSTART,IEND,JEND,X,Y,PHI,fid)

%---------------------- Capitulo 0 - Preliminares ------------------------%


%--------------- Capitulo 1 - Inicialização e Cebeçalhos -----------------%

ISKIP=1;
JSKIP=1;

%%%%% Alterável %%%%%
LINLIM=8; 

LINSTA=ISTART;

%Escreve o cabeçalho da matriz
fprintf(fid,'\r\n');
fprintf(fid,'--------------------------------------TEMPERATURE (ºC)---------------------------------------\r\n');

%Define um valor de LINEND inferior a IEND para inicializar o ciclo while
LINEND=0;

while LINEND<IEND
    
    %Escreve cabeçalho da matriz
    LINEND=LINSTA+(LINLIM-1)*ISKIP;
    LINEND=min(IEND,LINEND);
    
    %Escreve a segunda linha do output
    fprintf(fid,'I= ');
    for KI=LINSTA:ISKIP:LINEND
        fprintf(fid,'%10d',KI);
   
    end
    fprintf(fid,'       Y\r\n');
  
    
%----------------- Capitulo 2 - Escreve a matriz PHI ---------------------%
    
    fprintf(fid,'  J\r\n');
    for JJ=JSTART:JSKIP:JEND
        J=JSTART+JEND-JJ;
        IS=0;
        
        for I=LINSTA:ISKIP:LINEND
            A=PHI(I,J);
            
            IS=IS+1;
            
            if abs(A)<(1e-20)
                A=0.0;
            end
               
                STORE(IS)=A;
                
        end
        
        %Escreve o valor de J, T e Y
        fprintf(fid,'%3d',J);
        for KT=1:IS
        fprintf(fid,'%10.2E', STORE(KT));
        end
        fprintf(fid,'%10.4f\r\n', Y(J));
        
    end
        
    %Escreve uma linha com os X's
    fprintf(fid,'\r\n');
    fprintf(fid,'X= ');
    for KX=LINSTA:ISKIP:LINEND
        fprintf(fid,'%10.4f',X(KX));
    end
    
    fprintf(fid,'\r\n \r\n');
    
    LINSTA=LINEND+ISKIP;
   
end

fprintf(fid,'\r\n \r\n');

%Termina PRINT
end