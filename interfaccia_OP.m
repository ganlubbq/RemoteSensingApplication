%interfaccia dielettrica per onda piana (onda incidente e trasmessa)

clear all;
close all;
clc;

zlimu=30; %limite z up
zlimd=-60; %z limite down
xlimu=50;
xlimd=-50;
n=900;
[X,Z]=meshgrid(linspace(xlimu,xlimd,n),linspace(zlimu,zlimd,n)); %griglia di punti simmetrica

%strati
s1=+15; %localizzazione in z dell'inizio del secondo strato (circa, c'è il floor)
s2=-35; %del terzo strato (circa, c'è il floor)
%in pixel
spe1=floor((zlimu-s1)/((zlimu-zlimd)/n)); %spessore strato 1 in pixel
spe2=floor((zlimu-s2)/((zlimu-zlimd)/n)); %spessore strato 2 in pixel




%movie
%nframe=50; %number of time frames
t=0;
%for l=1:nframe


s=ones(size(X)); %matrice della velocità di propagazione del mezzo
e0=8.854*1e-12; %permettività elettrica del vuoto
mu0=4*pi*1e-7; %permeabilità magnetica del vuoto
c0=1/sqrt(e0*mu0); %velocità di propagazione della luce nel vuoto
eps=e0*ones(size(X));
mu=mu0*ones(size(X));
sigma=ones(size(X)); %matrice della conduttività dei mezzi (per mezzi reali)

%parametri dei mezzi

er1=6; %permettività dielettrica relativa, mezzo 1
er2=2; %permettività dielettrica relativa, mezzo 2
er3=3;
mur1=1; %perm magn relativa
mur2=1; %per materiali non magnetici è 1
mur3=1;
sigma1=0; %conduttività (zero per mezzi ideali)
sigma2=0;
sigma3=0;

c1=c0/sqrt(er1*mur1);
c2=c0/sqrt(er2*mur2);
c3=c0/sqrt(er3*mur3);

the=-10; %in gradi, angolo di incidenza, a partire dall'asse z, angoli positivi verso destra
the2=asin(sqrt((mur1*er1)/(mur2*er2))*sin(deg2rad(the))); %angolo di trasmissione dell' onda nel secondo mezzo
the3=asin(sqrt((mur2*er2)/(mur3*er3))*sin(the2));

%se c'è angolo critico, dimmelo

if(t==0)
thetac1=asin(sqrt(mod((mur2*er2)/(mur1*er1),pi)));
thetac2=asin(sqrt(mod((mur3*er3)/(mur2*er2),pi)));
disp('Angolo critico di ingresso nel mezzo 2 (gradi): '), rad2deg(thetac1)
disp('Angolo critico di ingresso nel mezzo 3 (gradi): '), rad2deg(thetac2)
end

%complichiamoci la vita
eps(1:spe1,:)=er1*eps(1:spe1,:);
eps(spe1+1:spe2,:)=er2*eps(spe1+1:spe2,:);
eps(spe2+1:end,:)=er3*eps(spe2+1:end,:);

mu(1:spe1,:)=mur1*mu(1:spe1,:);
mu(spe1+1:spe2,:)=mur2*mu(spe1+1:spe2,:);
mu(spe2+1:end,:)=mur3*mu(spe2+1:end,:);

sigma(1:spe1,:)=sigma1*sigma(1:spe1,:);
sigma(spe1+1:spe2,:)=sigma2*sigma(spe1+1:spe2,:);
sigma(spe2+1:end,:)=sigma3*sigma(spe2+1:end,:);




theta=ones(size(X)); %matrice angoli di incidenza
theta(1:spe1,:)=deg2rad(the)*ones(size(theta(1:spe1,:)));
theta(spe1+1:spe2,:)=the2*ones(size(theta(spe1+1:spe2,:)));
theta(spe2+1:end,:)=the3*ones(size(theta(spe2+1:end,:)));

f0=25*1e6; %frequenza temporale dell'onda monocromatica incidente, 25Mhz
dt=1/(3*f0); %per rispettare teo campionamento, basta 1/2*f0

w=2*pi*f0; %pulsazione temporale
k=w*sqrt(eps.*mu); %matrice del numero d'onda (nella direzione di propagazione)
kx=k.*sin(theta); %matrice della lentezza orizzontale
kz=sqrt(k.^2-kx.^2); %matrice della lentezza verticale

A=10; %ampiezza onda incidente
onda= A*exp(-j*kx.*X-j*kz.*Z-j*w*t);


%prova correzione di fase e ampiezza nel terzo strato (funziona)

%secondo strato
corz_f2= angle(onda(spe2,:))-angle(onda(spe2+1,:)); %correzione di fase da dare a tutti i valori dopo spe2
corz_a2=max(abs(onda(spe2,:)))/max(abs(onda(spe2+1,:)));
%primo strato
corz_f1= angle(onda(spe1,:))-angle(onda(spe1+1,:)); %correzione di fase da dare a tutti i valori dopo spe2
corz_a1=max(abs(onda(spe1,:)))/max(abs(onda(spe1+1,:)));

%primo strato
onda(spe1+1:spe2,:)=(corz_a1*abs(onda(spe1+1:spe2,:))).*(exp(j*(angle(onda(spe1+1:spe2,:))+repmat(corz_f1,spe2-spe1,1)))); %fase corretta

%secondo strato
onda(spe2+1:end,:)=(corz_a2*abs(onda(spe2+1:end,:))).*(exp(j*(angle(onda(spe2+1:end,:))+repmat(corz_f2,size(Z,1)-spe2,1)))); %fase corretta


%plotting


% h=mesh(Z,X,real(onda)),view([90 270]);
% colormap('gray');
% t=t+dt;
% set(h,'CData',real(onda))
% M(l)=getframe;
% end

figure,mesh(Z,X,real(onda)),view([90 270]),colormap('gray');


%play movie

%figure,movie(M,1,2)

