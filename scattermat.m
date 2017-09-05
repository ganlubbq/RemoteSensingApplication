%MAtrice di scatter per interfaccia elastica solido-solido

function [S,the2,thes,thes2,xi2,n2]=scattermat(a1,b1,ro1,a2,b2,ro2,thei)
%Parametri in ingresso (tutti parametri caratteristici del mezzo tranne l'angolo di incidenza)

% syms a1; %vel prop onde P mezzo 1
% syms a2; %vel prop onde P mezzo 2
% syms b1; %vel prop onde S mezzo 1
% syms b2; %vel prop onde S mezzo 2
% syms ro1; %densità mezzo 1
% syms ro2; %densità mezzo 2
% syms the; %angolo di incidenza onda P (va inserito il valore in radianti)
% 
% %Parametri derivati
% 
% syms thes; %angolo di incidenza onda S nel mezzo 1
% syms thes2; %angolo di incidenza onda S nel mezzo 2
% syms the2; %angolo di incidenza onda P nel mezzo 2
% syms p; %lentezza orizzontale, comune a tutte le onde, sin(the)/a1=sin(the2)/a2 e sin(thes)/b1=sin(thes2)/b2 in più sin(the)/a1=sin(thes)/b1
% % e sin(the2)/a2=sin(thes2)/b2
% syms n1; %lentezza verticale nel mezzo 1 onde S
% syms n2; %lentezza verticale nel mezzo 2 onde S
% syms xi1; %lentezza verticale nel mezzo 1 onde P
% syms xi2; %lentezza verticale nel mezzo 2 onde P
% syms w; %pulsazione temporale
% 
% %Onde Incidenti
% 
% syms Pd1;
% syms Sd1;
% syms Pu2;
% syms Su2;

%Onde scatterate
% syms Pu1;
% syms Su1;
% syms Pd2;
% syms Sd2;i

%thei = già in radianti

the2=asin((a2/a1)*sin(thei));
thes=asin((b1/a1)*sin(thei));
thes2=asin((b2/a2)*sin(the2));

p=sin(thei)/a1;

n1=cos(thes)/b1;

n2=cos(thes2)/b2;

xi1= cos(thei)/a1;

xi2= cos(the2)/a2;

%out-->M
M=[-xi1 p -xi2 -p; ...
    p n1 -p n2;...
    -ro1*(1-2*b1^2*p^2) 2*ro1*b1^2*p*n1 ro2*(1-2*b2^2*p^2) 2*ro2*b2^2*p*n2;
    2*ro1*b1^2*p*xi1 ro1*(1-2*b1^2*p^2) 2*ro2*b2^2*p*xi2 -ro2*(1-2*b2^2*p^2)];
%     -ro1*(b1^2/a1^2)*sin(2*thei) -ro1*(cos(thes)^2 - sin(thes)^2) -ro2*(b2^2/a2^2)*sin(2*the2) ro2*(cos(thes2)^2 - sin(thes2)^2);...%txz
%     ro1*(1-2*sin(thes)^2) -ro1*sin(2*thes) ro2*(-1+2*sin(thes2)^2) -ro2*sin(2*thes2)];%tzz
     
    
    
%IN --> N
N=[-xi1 -p -xi2 p;...
    p -n1 -p -n2;...
    ro1*(1-2*b1^2*p^2) 2*ro1*b1^2*p*n1 -ro2*(1-2*b2^2*p^2) 2*ro2*b2^2*p*n2;
    2*ro1*b1^2*p*xi1 -ro1*(1-2*b1^2*p^2) 2*ro2*b2^2*p*xi2 ro2*(1-2*b2^2*p^2)];
%     -ro1*(b1^2/a1^2)*sin(2*thei) ro1*(cos(thes)^2 - sin(thes)^2) -ro2*(b2^2/a2^2)*sin(2*the2) -ro2*(cos(thes2)^2 - sin(thes2)^2);...%txz
%     -ro1*(1-2*sin(thes)^2) -ro1*sin(2*thes) -ro2*(-1+2*sin(thes2)^2) -ro2*sin(2*thes2)];%tzz


%angoli di incidenza nel primo e secondo mezzo per onde P ed onde S


the2=subs(thei);
thes=subs(thei);
thes2=subs(the2);

%calcola la matrice di scatter
M=subs(M);
N=subs(N);

S=inv(M)*N;
