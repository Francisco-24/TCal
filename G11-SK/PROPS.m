%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                      script: PROPS (vers�o MATLAB)
%  Especifica��o da varia��o das propriedades do material, nomeadamente a
%                        condutividade t�rmica, K
%                                FS
%                            09/05/2011
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

function[GAMH]=PROPS(NI, NJ, TCON, GAMH)

%-----------------------Capitulo 0 - Preliminares-------------------------%


%------------- Capitulo 1 - Actualiza��o das propriedades ----------------%

for I=1:NI;
    for J=1:NJ;
        GAMH(I,J)=TCON(I,J);
    end
end


%Termina PROPS
end