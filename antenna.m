%diagramma di radiazione di un apertura rettangolare

clear all;
close all;
clc;

f0= 50*1e6; %frequenza operativa dell'antenna, 250 Mhz
c=3*1e8; %m/s velocità di propagazione della luce nel vuoto
lam= c/f0; %lunghezza d'onda
a= 20; %dimensione in metri della base dell'apertura, in x (verticale)
b= 10; %dimensione in metri della base dell'apertura in y (orizzontale)

%angoli per i quali viene calcolato il diagramma di radiazione, a theta
%fisso

%the=linspace(-pi/2,pi/2,400);
T=linspace(-pi/2,pi/2,2000);
ph=deg2rad(0); %campionamento angolare
%[T P]=meshgrid(the,ph); %griglia di campionamento angolare
fn= sinc(a*((sin(T).*cos(ph))./lam)).^2 .* sinc(b*((sin(T).*sin(ph))./lam)).^2;

figure,plot(T,10*log10(fn)),title(strcat('diagramma di radiazione per phi= ',num2str(rad2deg(ph)))),xlabel('Theta'),ylabel('Guadagno in dB')

