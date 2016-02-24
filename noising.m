function [signal_noised,variance_bruit,Pbwi]=noising(noise,signal,RSB)

fe=8000;
fo=3000;
N=length(signal);
t = (1:N)/fe ;

Bi=transpose(randn(1,N));
Sin=transpose(sin(2*pi*fo*t));

Ps = mean(abs(signal).^2); 
Pbwi = mean(abs(Bi).^2); 
Pbc=mean(abs(Sin).^2);

variance_bruit=(Ps/Pbwi)*10^(-RSB/10);

alpha_c=sqrt((Ps/Pbc)*10.^(-RSB/10));
B=sqrt(variance_bruit)*Bi;

if(noise=='w')
    
    signal_noised=signal+B;
    
else if(noise=='c')
        
        signal_noised=signal+alpha_c*Sin;
        
    end 
end
    

    
