%prova raggi geometrici

clear all;
close all;
clc;

x=linspace(0,6000,1000);
z=linspace(0,3000,1000);




dx=(max(x)-min(x))/length(x);
dz=(max(z)-min(z))/length(z);

%modello di velocità

c0=4000;
C=ones(length(z),length(x))*c0-repmat(0.5*z',1,length(x));

S=1./C; %lentezza

%posizione sorgente

x0=length(x)/2;
z0=0;


figure, imagesc(C),colormap(gray),colorbar, hold on, plot(z0,x0,'*r'),

%tracciamento
i=1;
for th=-13:1:13
i=1;
clear the;
the(i)=deg2rad(th); %angolo di start

xp(i)=x0;
zp(i)=z0;
px=sin(the(1))*S(i,1); %s/m in x 
pz(i)=cos(the(i)); %s/m in z
i=2;
while(i<=length(z))
the(i)=real(asin(px./S(i,1)));
pz(i)=cos(the(i))*S(i,1);

xp(i)=xp(i-1)+sin(the(i))*sqrt(dx^2+dz^2);
zp(i)=zp(i-1)+cos(the(i))*sqrt(dx^2+dz^2);
i=i+1;
end


%plotting
hold on, line(xp,zp), 
end