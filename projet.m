clear all;
close all;

fe=8000;
load('fcno01fz.mat');
load('fcno02fz.mat');
load('fcno03fz.mat');
load('fcno04fz.mat');
load('fcno05fz.mat');
N=512;

% Représentation temporelle et spectrogramme du signal original
figure, 
subplot(2,1,1)
plot(fcno01fz)
axis([0 length(fcno01fz) min(fcno01fz) max(fcno01fz)])
xlabel('temps')
ylabel('Amplitude')
title('Représentation temporelle du signal original')

subplot(2,1,2)
spectrogram(fcno01fz,N,0,N,fe,'yaxis')
title('Spectrogramme du signal original')

% Représentation temporelle et spectrogramme du signal bruité ('w' : bruité
% par un bruit blanc , 'c' : bruité par une sinusoïde)
[signal_noised,variance_bruit]=noising('w',fcno01fz,15);

figure, 
subplot(2,1,1)
plot(signal_noised)
axis([0 length(signal_noised) min(signal_noised) max(signal_noised)])
xlabel('temps')
ylabel('Amplitude')
title('Représentation temporelle du signal bruité')

subplot(2,1,2)
spectrogram(signal_noised,N,0,N,fe,'yaxis')
title('Spectrogramme du signal bruité')

%% Débruitage du signal bruité

%% 1ère méthode : fft du signal bruité en totalité, mise en place des 
%composantes fréquentielles associées au bruit à 0 puis transformée de 
%Fourier inverse de la séquence résultante
TFD_noise=fft(signal_noised);   %Calcul de la transformée de Fourier 
%de la totalité du signal de parole bruité

%On met les composantes fréquentielles associées au bruit à 0
TFDmax=max(abs(TFD_noise));        
[xmax]=find(abs(TFD_noise)==TFDmax);
TFD_noise(xmax)=0;

% On prend la transformée de Fourier inverse de la séquence résultante
%Puis on trace la réprésentation temporelle et le spectrogramme du signal
%débruité pour les comparer avec celles du signal original
signal_unnoised=ifft(TFD_noise);
figure,
subplot(2,1,1)
plot(signal_unnoised)
axis([0 length(fcno01fz) min(fcno01fz) max(fcno01fz)])
xlabel('temps')
ylabel('Amplitude')
title('Représentation temporelle du signal débruité (Méthode 1)')

subplot(2,1,2)
spectrogram(signal_unnoised,N,0,N,fe,'yaxis')
title('Spectrogramme du signal débruité (Méthode 1)')

%% 2ème méthode : Approche avec un filtre RIF/RII

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
title('Représentation temporelle du signal débruité (Méthode 2 RIF)')
xlabel('Temps')
ylabel('Amplitude')

subplot(2,1,2)
spectrogram(signal_unnoised_RIF,N,0,N,fe,'yaxis');
title('Spectrogramme du signal débruité (Méthode 2 RIF)');

%Cas du RII
rz=1; % Rayon des zéros
z1=rz*exp(-j*2*pi*f0/fe); % f0 fréq à supprimer & fe fréq échantillonnage
z2=conj(z1);
rp=0.9; % Rayon des pôles
p1=rp*exp(-j*2*pi*f0/fe); 
p2=conj(p1);
num=[1 -(z1+z2) rz^2];
denom=[1 -(p1+p2) rp^2];
signal_unnoised2=filter(num,denom,signal_noised);

figure,

subplot(2,1,1)
plot(signal_unnoised2)
axis([0 length(signal_unnoised2) min(signal_unnoised2) max(signal_unnoised2)])
title('Représentation temporelle du signal débruité (Méthode 2 RII)')
xlabel('Temps')
ylabel('Amplitude')

subplot(2,1,2)
spectrogram(signal_unnoised2,N,0,N,fe,'yaxis');
title('Spectrogramme du signal débruité (Méthode 2 RII)');

%% Procédure d'addition-recouvrement + Rehaussement : Méthode dite en sous-espaces

[signal_reconstructed,RSB_final]=Overlap_Add(signal_noised,variance_bruit);

figure, 
subplot(2,1,1)
plot(signal_reconstructed)
axis([0 length(signal_reconstructed) min(signal_reconstructed) max(signal_reconstructed)])
xlabel('temps')
ylabel('Amplitude')
title('Représentation temporelle du signal rehaussé')

subplot(2,1,2)
spectrogram(signal_reconstructed,N,0,N,fe,'yaxis')
title('Spectrogramme du signal rehaussé')





