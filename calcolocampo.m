%Calcolo del campo
clear all;
close all;
clc;



RdPP1 = struct2array(load('RdPP1_200.mat'));  
RdPS1 = struct2array(load('RdPS1_200.mat')); %solo riflessioni

RdPP2 = struct2array(load('RdPP2_200.mat'));
RdPS2 = struct2array(load('RdPS2_200.mat'));
TdPP1 = struct2array(load('TdPP1_200.mat'));


dthe=88.75/356; %risoluzione in angolo


dfR=50/1000; %risoluzione in frequenza della riflettività matrice


%forma d'onda in tx: Ricker : ho cambiato df=1Hz

fc=100;
dt=1/fc; %intervallo cmapionamento temporale
dur=10; %durata in secondi
t=-dur/2:dt:dur/2;
np=round(dur/dt); %numero di punti totale
if(mod(np,2)==0)
    np=np+1;
end
lat=(np-1)/2;
ts=0; %time shift
cfr=5; %frequenza centro banda

s=10*ricker(np,cfr,dt,ts,lat);

%Nfft=1000;
Nfft=length(s);
df=fc/Nfft; %risoluzione in frequenza della fft della Ricker
%t=-dur/2:dt:dur/2;
fk=assefr(Nfft,fc);
S=fft(s,Nfft);
%San=[1/2*S(502) S(503:1002)]; %501 campioni

%estraggo segnale analitico della forma d'onda
San=[1/2*S(1) S(2:(Nfft-1)/2+1)];

%Parametri mezzo 1

a1=2000; %m/s mezzo 1 onda P
b1=1329; %m/s mezzo 1 onda S
a2=3000;
b2=1754;

%distanze geofoni
rr=50:25:2500; %11 rx
zr=0;

h1=600; %distanza sorgente - interfaccia
h2=200;

l1=sqrt(4*h1^2+rr.^2); %distanze dei raggi
l2=sqrt(4*h2^2+rr.^2) + l1;
%distanze raggio di Fermat (ritardi in campioni sorgente-ricevitore)
drfP1=((sqrt((2*h1)^2+(rr).^2))./a1);
drfS1=(sqrt(h1^2+(rr/2).^2)./a1+sqrt(h1^2+(rr/2).^2)./b1);

drfP2=((sqrt((2*h2)^2+(rr).^2))./a2) + drfP1;
drfS2=(sqrt(h2^2+(rr/2).^2)./a2+sqrt(h2^2+(rr/2).^2)./b2) +drfS1;


%   Integrazione con discretizzazione angolare e in frequenza
%
% %f=0.1:0.1:50;
% f=0:df:50; %campioni in frequenza per i quali abbiamo la funz di riflettività
% the=0.25:0.25:38; %superiamo angolo critico, mettendo piccola parte immaginaria a kz, partiva da 0.25
% 
% for q=1:length(rr)
%     clear phu;
%     clear psu;
% for v=1:length(f) %cicla su tutte le frequenze 
%     
%     w=2*pi*f(v);
%     Jf=0;
%     Yf=0;
%     for m=1:length(the) %cicla su tutte le direzioni, per ogni frequenza
%     
%     
%     thei=deg2rad(the(m)); %per questo angolo
% 
% 
%     %kzp=w*cos(thei)/a + + j*((w*cos(thei)/a)/1e12);
%     %kzs=w*cos(thei)/b + j*((w*cos(thei)/b)/1e12);
%     %krp=sqrt((w/a)^2-kzp^2);
%     %krs=sqrt((w/b)^2-kzs^2);
%     p=sin(thei)/a; %lentezza radiale
%     krp=w*p +  1e-16;
%      %j*((w*sin(thei)/b)/1e12);
%     xi=sqrt((1/a)^2-p^2);
%     kzp=sqrt((w/a)^2-krp^2);
%     kzs=sqrt((w/b)^2-krp^2);
%     
%     %valutazione dell'integrale con forza bruta
%     
%     Jf=Jf + (besselj(0,krp*rr(q))/kzp).*RdPP1(floor(rad2deg(thei)/dthe),floor(((w/(2*pi))/df)+1))*exp(-j*kzp*(2*h-zr)).*krp;
%     Yf=Yf + (besselj(0,krp*rr(q))/kzs).*RdPS1(floor(rad2deg(thei)/dthe),floor(((w/(2*pi))/df)+1))*exp(-j*(kzp*h+kzs*(h-zr))).*krp;
%        
%     end
%     
%     %Jf=[conj(Jf(501:-1:2)) Jf];
%     %Yf=[conj(Yf(501:-1:2)) Yf];
%     Jf=j*Jf;
%     Yf=j*Yf;
%     
%     phu(v)=San(v)*Jf;
%     psu(v)=San(v)*Yf;
% 
% end
% % 
% % %integrali completi, per ogni angolo --> riga, dunque per ogni onda conica,
% % %ho la risposta in frequenza della stratificazione, per ottenere la
% % %completa devo sommare tutte le onde
% % 
% ritP=exp(-j*2*pi*[0:(Nfft-1)/2]*df*(drfP(q)));
% ritS=exp(-j*2*pi*[0:(Nfft-1)/2]*df*(drfS(q)));
% 
% phu=phu.*ritP;
% psu=psu.*ritS;
% 
% phu=[phu conj(phu(2:(Nfft-1)/2+1))];
% psu=[psu conj(psu(2:(Nfft-1)/2+1))];
% 
% sisP(:,q)=ifft(phu,'symmetric');
% sisS(:,q)=ifft(psu,'symmetric');
% end

%valutazione integrale con il metodo della fase stazionaria

%non ciclo più sugli angoli, poichè per ogni ricevitore guardo solo le
%direzioni che contribuiscono in fase rispetto al raggio di Fermat, che è
%il raggio corrispondente al parametro p0, punto a fase stazionaria

%Integrazione con il metodo della fase stazionaria

f=0.1:df:50.1;

for q=1:length(rr) %per ogni ricevitore
clear temp;
clear phu;
for v=1:length(f)
    w=2*pi*f(v);
    %onde P nel mezzo 1
    p0=rr(q)/(a1*l1(q));
    thei=asin(p0*a1);
    xi0=(2*h1)/(a1*l1(q));
    phu1(v)=(sqrt(2/(pi*w*p0*rr(q)))*(RdPP1(floor(rad2deg(thei)/dthe),floor(((w/(2*pi))/df)+1))/(j*xi0))*(w*p0)*(sqrt(2*pi/(w*((a1*l1(q)^3)/4*h1^2))))*exp(-j*w*(l1(q)/a1)+(pi/4)));
    
    %onde S nel mezzo 1
    p0=rr(q)/(b1*l1(q));
    thei=asin(p0*b1);
    ni0=(2*h1)/(b1*l1(q));
    psu1(v)=(sqrt(2/(pi*w*p0*rr(q)))*(RdPS1(floor(rad2deg(thei)/dthe),floor(((w/(2*pi))/df)+1))/(j*ni0))*(w*p0)*(sqrt(2*pi/(w*((b1*l1(q)^3)/4*h1^2))))*exp(-j*w*(l1(q)/b1)+(pi/4)));
    
    %ricevitori nello strato 2, onde P nel mezzo 2
    p0=rr(q)/(a2*l2(q));
    thei=asin(p0*a2);
    xi0=(2*h2)/(a2*l2(q));
    phu2(v)=(sqrt(2/(pi*w*p0*rr(q)))*(RdPP2(floor(rad2deg(thei)/dthe),floor(((w/(2*pi))/df)+1))/(j*xi0))*(w*p0)*(sqrt(2*pi/(w*((a2*l2(q)^3)/4*h2^2))))*exp(-j*w*(l2(q)/a2)+(pi/4)));
    
    %onde S nel mezzo 2
    p0=rr(q)/(b2*l2(q));
    thei=asin(p0*b2);
    ni0=(2*h2)/(b1*l2(q));
    psu2(v)=(sqrt(2/(pi*w*p0*rr(q)))*(RdPS2(floor(rad2deg(thei)/dthe),floor(((w/(2*pi))/df)+1))/(j*ni0))*(w*p0)*(sqrt(2*pi/(w*((b2*l2(q)^3)/4*h2^2))))*exp(-j*w*(l2(q)/b2)+(pi/4)));
    
end

%onde P nel mezzo 1
temp1=San.*phu1;
rit1=exp(-j*2*pi*[0:length(f)-1]*df*drfP1(q));
temp1=temp1.*rit1;
temp1=[temp1 conj(temp1(2:(Nfft-1)/2+1))];
sisPstaz1(:,q)=ifft(temp1,'symmetric');


%onde S nel mezzo 1
temp1s=San.*psu1;
rit1s=exp(-j*2*pi*[0:length(f)-1]*df*drfS1(q));
temp1s=temp1s.*rit1s;
temp1s=[temp1s conj(temp1s(2:(Nfft-1)/2+1))];
sisSstaz1(:,q)=ifft(temp1s,'symmetric');

%nel mezzo 2 entra la forma d'onda propagata fino allo strato1 e trasmessa,
%onda P

temp2=San.*phu2;
rit2=exp(-j*2*pi*[0:length(f)-1]*df*drfP2(q));
temp2=temp2.*rit2;
temp2=[temp2 conj(temp2(2:(Nfft-1)/2+1))];
sisPstaz2(:,q)=ifft(temp2,'symmetric');

%onde S nel mezzo 2
temp2s=San.*psu2;
rit2s=exp(-j*2*pi*[0:length(f)-1]*df*drfS2(q));
temp2s=temp2s.*rit2s;
temp2s=[temp2s conj(temp2s(2:(Nfft-1)/2+1))];
sisSstaz2(:,q)=ifft(temp2s,'symmetric');

end


%Plotting

seismic=s_convert(sisPstaz1,0,dur);
s_wplot(seismic,{'scale','no'},{'peak_fill',''})


seismic2=s_convert(sisPstaz2,0,dur);
s_wplot(seismic2,{'scale','no'},{'peak_fill',''})

seismicS=s_convert(sisSstaz1,0,dur);
s_wplot(seismicS,{'scale','no'},{'peak_fill',''})


seismic2S=s_convert(sisSstaz2,0,dur);
s_wplot(seismic2S,{'scale','no'},{'peak_fill',''})