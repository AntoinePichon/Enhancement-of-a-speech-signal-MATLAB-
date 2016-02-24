function [signal_reconstructed,RSB_final]=Overlap_Add(signal,variance_bruit)

fe=8000;
t_frame=0.025;
N_ech=t_frame*fe;

    signal(length(signal)+1:length(signal)+((round(length(signal)/N_ech)+1)*N_ech)-length(signal))=0;
    
    signal_pad=signal;
    
    %signal_framed=reshape(signal_pad,N_ech,length(signal_pad)/N_ech);
    
covering=50/100;
length_covering=covering*N_ech;
Nb_frame=length(signal)/length_covering-1;

window=hamming(N_ech);

signal_framed=zeros(Nb_frame,N_ech);
concat_frame=zeros(1,length(signal));
concat_window=zeros(1,length(signal));

%% Decoupage du signal en trames

        for i=1:Nb_frame
            signal_framed(i,:)=signal_pad(((i-1)/2)*N_ech+1:((i-1)/2)*N_ech+N_ech);
            
            signal_framed(i,:)=signal_framed(i,:).*transpose(window);
        end
        
       %% Rehaussement 
       signal_framed=raising(signal_framed,variance_bruit);
       
        %% Reconstruction
        
        for k=1:Nb_frame
            
            concat_frame(((k-1)/2)*N_ech+1:(((k-1)/2)+1)*N_ech)= concat_frame(((k-1)/2)*N_ech+1:(((k-1)/2)+1)*N_ech)+signal_framed(k,:);
            
            concat_window(((k-1)/2)*N_ech+1:(((k-1)/2)+1)*N_ech)= concat_window(((k-1)/2)*N_ech+1:(((k-1)/2)+1)*N_ech)+transpose(window);
            
        end
        
        signal_reconstructed=concat_frame./concat_window;
            
        signal_reconstructed(isnan(signal_reconstructed))=0;
        
      Ps_final = mean(abs(signal_reconstructed).^2);

      Pb_final = mean(abs(signal_pad-transpose(signal_reconstructed)).^2);

      RSB_final=10*log10(Ps_final/Pb_final);
        
    end
    
    





