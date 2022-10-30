%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: PROPS (versão MATLAB)
%  Especificação da variação das propriedades do material, nomeadamente a
%                        condutividade térmica, K
%                                FS
%                            09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[GAMH]=PROPS(NI, NJ, TCON, GAMH)

%-----------------------Capitulo 0 - Preliminares-------------------------%


%------------- Capitulo 1 - Actualização das propriedades ----------------%

for I=1:NI;
    for J=1:NJ;
        GAMH(I,J)=TCON(I,J);
    end
end


%Termina PROPS
end