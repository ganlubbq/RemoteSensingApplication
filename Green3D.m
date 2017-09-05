%Prova funzioni di Green 3D

clear all;
clc;
close all;

x=linspace(0,3,101);
y=linspace(0,8,101);
z=linspace(0,8,101);

[X,Y,Z]=meshgrid(x,y,z);

% R=sqrt(X.^2+Y.^2+Z.^2); %matrice che per ogni punto, contiene il raggio della sfera corrispondente
% 
% %rappresentazione in (r,w)
% 
 c=300; %velocità del mezzo
% w=30; %rad/s
% 
% P=1./(4*pi*R).*exp(j*w*(R/c)); %funz di Green 3D
% 
% %Visione
% 
% figure,imagesc(angle(P(:,:,1))),title('Fase')
% figure,imagesc(abs(P(:,:,1))),title('Ampiezza')

%Sorgente diversa da impulso, monocromatica a freq.30Hz e presente come 2 impulsi

% %coord 1 punto
% x01=0;
% y01=3.98;
% z01=2;
% f1=300;
% w1=2*pi*f1;
% 
% R01=sqrt((X-x01).^2+(Y-y01).^2+(Z-z01).^2);
% 
% P1=1./(4*pi*R01).*exp(j*w1*(R01/c)); %ho integrato su tutto il tempo, campo generato in tutto lo spazio, per 1 frequenza, dal punto1
% 
% %coord. punto 2
% x02=0;
% y02=6.02;
% z02=2;
f2=3000;
w2=2*pi*f2;

% R02=sqrt((X-x02).^2+(Y-y02).^2+(Z-z02).^2);
% 
% P2=1./(4*pi*R02).*exp(j*w2*(R02/c));

%Posso mettere la creazione dei punti in un ciclo, e generare il campo
%diffratto da un'apertura rettangolare


x0=0;
P=zeros(101,101,101);

for y0=4:.2:5
    for z0=4:.2:5

        R0=sqrt((X-x0).^2+(Y-y0).^2+(Z-z0).^2);
        P0=1./(4*pi*R0).*exp(j*w2*(R0/c));
        P=P+P0;
    end
end

%campo totale come antitrasformata (lo vedo in (r,t))
t=0;
%dt=0.005;
%for k=1:30
P=P0*exp(-j*w2*t);
imagesc(abs(P(:,:,2)).^2),colormap('gray')
%M(k)=getframe;
%t=t+dt;
%end
%close all;
%movie(M,1,2)