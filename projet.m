clear all;
close all;

fe=8000;
load('fcno01fz.mat');
load('fcno02fz.mat');
load('fcno03fz.mat');
load('fcno04fz.mat');
load('fcno05fz.mat');
N=512;

% Repr�sentation temporelle et spectrogramme du signal original
figure, 
subplot(2,1,1)
plot(fcno01fz)
axis([0 length(fcno01fz) min(fcno01fz) max(fcno01fz)])
xlabel('temps')
ylabel('Amplitude')
title('Repr�sentation temporelle du signal original')

subplot(2,1,2)
spectrogram(fcno01fz,N,0,N,fe,'yaxis')
title('Spectrogramme du signal original')

% Repr�sentation temporelle et spectrogramme du signal bruit� ('w' : bruit�
% par un bruit blanc , 'c' : bruit� par une sinuso�de)
[signal_noised,variance_bruit]=noising('w',fcno01fz,15);

figure, 
subplot(2,1,1)
plot(signal_noised)
axis([0 length(signal_noised) min(signal_noised) max(signal_noised)])
xlabel('temps')
ylabel('Amplitude')
title('Repr�sentation temporelle du signal bruit�')

subplot(2,1,2)
spectrogram(signal_noised,N,0,N,fe,'yaxis')
title('Spectrogramme du signal bruit�')

%% D�bruitage du signal bruit�

%% 1�re m�thode : fft du signal bruit� en totalit�, mise en place des 
%composantes fr�quentielles associ�es au bruit � 0 puis transform�e de 
%Fourier inverse de la s�quence r�sultante
TFD_noise=fft(signal_noised);   %Calcul de la transform�e de Fourier 
%de la totalit� du signal de parole bruit�

%On met les composantes fr�quentielles associ�es au bruit � 0
TFDmax=max(abs(TFD_noise));        
[xmax]=find(abs(TFD_noise)==TFDmax);
TFD_noise(xmax)=0;

% On prend la transform�e de Fourier inverse de la s�quence r�sultante
%Puis on trace la r�pr�sentation temporelle et le spectrogramme du signal
%d�bruit� pour les comparer avec celles du signal original
signal_unnoised=ifft(TFD_noise);
figure,
subplot(2,1,1)
plot(signal_unnoised)
axis([0 length(fcno01fz) min(fcno01fz) max(fcno01fz)])
xlabel('temps')
ylabel('Amplitude')
title('Repr�sentation temporelle du signal d�bruit� (M�thode 1)')

subplot(2,1,2)
spectrogram(signal_unnoised,N,0,N,fe,'yaxis')
title('Spectrogramme du signal d�bruit� (M�thode 1)')

%% 2�me m�thode : Approche avec un filtre RIF/RII

f0=3000;

% Filtrage du bruit

%Cas du RIF
r=1;
z=r*exp(-j*2*pi*f0/fe);
zconj=conj(z);
RIF=[1 -(z+zconj) r^2];

signal_unnoised_RIF=filter(RIF,1,signal_noised);

figure,

subplot(2,1,1)
plot(signal_unnoised_RIF)
axis([0 length(signal_unnoised_RIF) min(signal_unnoised_RIF) max(signal_unnoised_RIF)])
title('Repr�sentation temporelle du signal d�bruit� (M�thode 2 RIF)')
xlabel('Temps')
ylabel('Amplitude')

subplot(2,1,2)
spectrogram(signal_unnoised_RIF,N,0,N,fe,'yaxis');
title('Spectrogramme du signal d�bruit� (M�thode 2 RIF)');

%Cas du RII
rz=1; % Rayon des z�ros
z1=rz*exp(-j*2*pi*f0/fe); % f0 fr�q � supprimer & fe fr�q �chantillonnage
z2=conj(z1);
rp=0.9; % Rayon des p�les
p1=rp*exp(-j*2*pi*f0/fe); 
p2=conj(p1);
num=[1 -(z1+z2) rz^2];
denom=[1 -(p1+p2) rp^2];
signal_unnoised2=filter(num,denom,signal_noised);

figure,

subplot(2,1,1)
plot(signal_unnoised2)
axis([0 length(signal_unnoised2) min(signal_unnoised2) max(signal_unnoised2)])
title('Repr�sentation temporelle du signal d�bruit� (M�thode 2 RII)')
xlabel('Temps')
ylabel('Amplitude')

subplot(2,1,2)
spectrogram(signal_unnoised2,N,0,N,fe,'yaxis');
title('Spectrogramme du signal d�bruit� (M�thode 2 RII)');

%% Proc�dure d'addition-recouvrement + Rehaussement : M�thode dite en sous-espaces

[signal_reconstructed,RSB_final]=Overlap_Add(signal_noised,variance_bruit);

figure, 
subplot(2,1,1)
plot(signal_reconstructed)
axis([0 length(signal_reconstructed) min(signal_reconstructed) max(signal_reconstructed)])
xlabel('temps')
ylabel('Amplitude')
title('Repr�sentation temporelle du signal rehauss�')

subplot(2,1,2)
spectrogram(signal_reconstructed,N,0,N,fe,'yaxis')
title('Spectrogramme du signal rehauss�')





