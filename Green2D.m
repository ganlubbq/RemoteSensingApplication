%Prova funz di Green2D, diffrazione da un'apertura SLIT

clear all;
close all;
clc;

z=linspace(0,20,1024); %direz verticale, dove ho l'apertura
dz=(max(z)-min(z))/length(z);
y=linspace(0,60,1024); %direz orizzontale
dy=(max(y)-min(y))/length(y);

fcz=1/dz;
fcy=1/dy;

[Y,Z]=meshgrid(y,z);

c=300;
f=3000; %dunque lambda 1

lam=c/f;

w=2*pi*f;

y0=0;

P1=zeros(size(Y));

zin=10;

%P2=zeros(size(Y));

for z0=zin:lam/5:zin+10*lam  %apertura slit da 5*lambda

    R0=sqrt((Y-y0).^2+(Z-z0).^2);
    P0=5*j/(4*pi)*besselh(0,1,w*R0/c); %funz Green 2d
    %P02=sqrt(2./(pi*R0)).*exp(R0-(c/2)*pi-pi/4);
    P1=P1+P0;
    %P2=P2+P02;
end

figure,imagesc(y,z,(real(P1).^2+imag(P1).^2)),colormap(gray),caxis([0 0.5]),title('Energia del campo ')

Fp1=fftshift(fft2(P1));
fky=0:fcy/length(y):fcy/2-fcy/length(y);
fkz=-fcz/2:fcz/length(z):fcz/2-fcz/length(z);

figure,imagesc(fky,fkz,abs(Fp1(:,1:length(fky)))),colormap(gray),title('Spettro nel dominio dei numeri d''onda')
figure,imagesc(y,z,real(P1)),colormap(gray),title('Campo')
%figure,imagesc(real(P2).^2+imag(P2).^2),colormap('gray'),title('Hankel approx')