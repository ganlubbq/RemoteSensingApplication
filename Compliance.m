% Compliance matrix e stiffness matrix col Symbolic Matlab
clear all;
close all;
clc;

syms E; %modulo di Young
syms v; %coeff di Poisson
S=sym(1/E.*[1 -v -v 0 0 0;-v 1 -v 0 0 0;-v -v 1 0 0 0;0 0 0 1+v 0 0;0 0 0 0 1+v 0;0 0 0 0 0 1+v]); %stiffness isotropa
C=inv(S);

%sostituzione di valori numerici nelle matrici simboliche così calcolate
y=2*1e11; %modulo di Young di un materiale, acciaio
p=0.33; %modulo di Poisson di un materiale, acciaio
Cnum=subs(C,{E,v},{y,p});
Snum=subs(S,{E,v},{y,p});
disp('Matrice di Compliance per acciaio isotropo')
Cnum
disp('Matrice di Stiffness per acciaio isotropo')
Snum

sf=zeros(6,1); %tensore sforzi in notazione abbreviata
def=zeros(6,1); %tensore deformazioni in notazione abbreviata

%Imprimo uno sforzo normale in direzione x
%imprimo sforzo normale gradualmente crescente, fino a 1e5 [N/m^2]
i=1;
sigx(i)=10;
sf(1)=sigx(i); %[N/m^2]
while (sigx(i)<=1e5)
def(:,i)=Snum*sf;
sf(1)=sf(1)+sigx(i); %aumento di 10 Newton ogni volta
sigx(i+1)=sigx(i)+10;
i=i+1;
end
sumdef=zeros(1,length(def(1,:)));
sumdef=def(1,:)+def(2,:)+def(3,:); %mostro come varia la deformazione volumetrica dell'acciaio, all'aumentare dello sforzo normale
figure,plot(sigx(1:length(sumdef)),sumdef,'k'),title('Andamento della deformazione volumetrica dell''acciaio in funzione di uno sforzo normale applicato in direzione x'),...
    xlabel('sigmax [N/m^2]'),ylabel('Deformazione'),...
    hold on,...
    plot(sigx(1:length(sumdef)),def(1,:),'r-.'),...
    plot(sigx(1:length(sumdef)),def(2,:),'b-.'),legend('Deformazione volumetrica','Deformazione assiale','Deformazione nel piano normale')

