%Riflettività 1D

clear all;
close all;
clc;

%Parametri dei mezzi e operativi

% f=10;
% w=2*pi*f;

a1=2000;
b1=1329;
ro1=2;
a2=3000;
b2=1754;
ro2=2.7;
a3=2000;
b3=1332;
ro3=2;



%for z=1:3
%thei=30;
h2=50; %spessore strato 2

i=1;
for thei=90:-0.25:0
    k=1;
    the=deg2rad(thei);
    [S1,the2,thes,thes2,xi2,n2]=scattermat(a1,b1,ro1,a2,b2,ro2,the); %matrice di scatter per interfaccia 1 ed angoli di :

        R1D=S1(1:2,1:2);
        T1U=S1(1:2,3:4);
        T1D=S1(3:4,1:2);
        R1U=S1(3:4,3:4);

% thes --> incidenza S
% the2 --> tx P
% thes2 --> tx S


     [S2,the3,thes2,thes3,xi3,n3]=scattermat(a2,b2,ro2,a3,b3,ro3,the2); %matrice di scatter interf 2

        R2D=S2(1:2,1:2);
        T2U=S2(1:2,3:4);
        T2D=S2(3:4,1:2);
        R2U=S2(3:4,3:4);

    
    for f=0:0.1:50;
        w=2*pi*f;
        
%onde incidenti dall'alto

        O1d=[1;1];
        Q2=[exp(j*w*xi2*h2) 0;0 exp(j*w*n2*h2)]; %propagatore per il mezzo 2

%funz riflettività per interf 1
        RR2D=R2D;
        RR1D=R1D+T1U*Q2*RR2D*Q2*(diag(size(Q2))+(R1U*Q2*RR2D*Q2))*T1D; %1 riverberazione

%calcolo campo up nel mezzo 1
%         O1u=RR1D*O1d;
        RifPPd(i,k)=RR1D(1,1);
        
        k=k+1;
    end
    i=i+1;
end

figure,imagesc(abs(RifPPd)),colormap(gray)
