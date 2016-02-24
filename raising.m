function [signal_raised]=raising(signal_framed,variance_bruit)

fe=8000;
t_frame=0.025;
N_ech=t_frame*fe;

[Nb_frame Size_frame]=size(signal_framed);

L=70;

M=Size_frame+1-L;

for i=1:Nb_frame
    
    one_frame=signal_framed(i,:);
    
    H=hankel(one_frame(1:L),one_frame(L:Size_frame));
    
    [U,S,V]=svd(H);
    
    seuil=sqrt(variance_bruit*M);
    
    [row, column]=size(S);
    
    k=1;
    while k<=column && S(k,k)>=seuil;
        k=k+1;
    end
    
    K=k-1;
    
    S_dom=S(1:K,1:K);
    
    U_dom=U(:,1:K);
    V_dom=V(:,1:K);
    
    H_dom=U_dom*S_dom*V_dom';
    
    for j=-(L-1):M-1
        
        tmp(j+L)=mean(diag(fliplr(H_dom),j));
    end
    
    signal_raised(i,:)=fliplr(tmp);
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    