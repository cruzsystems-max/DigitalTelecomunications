function [gray,bk,fik,Ik,Qk]= table_cor(n,N)%Tabla de correspondecia
    M=2^n; 
    bk=0:1:2^n-1;
    dec=decTobin(bk);
    gray=bin2gray(dec);
    fik=pi*(2*bk+N)/M; 
    Ik=cos(fik);
    Qk=sin(fik);
end

