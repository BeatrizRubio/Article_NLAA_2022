%Cálculo de la factorización de la matriz wronkiana de la base de
%Bernstein con la factorizacion LDU (resulta ser la misma que BDA con los multipliers)

clear all

format long E
 syms t
 
n=4

%  syms x
 
 f=t
 g=1-t
 f1=cos(2*pi*t)
 f2=sin(2*pi*t)
 t=0.5

 %Calculo de las x


for i=1:n+1
	xb(i)=eval( diff(f1,i-1)); 	
end
xb
for i=1:n+1
	yb(i)=eval( diff(f2,i-1)); 	
end

% %Cálculo de la matriz wronskiana
% 
 A=zeros(n+1);

for i=1:n+1
	for j=1:n+1
	A(i,j)=eval( diff(nchoosek(n,j-1)*  ((f))^(j-1) *((g))^(n+1-j),i-1));
 	end 
end





% 
% %Cálculo de los multiplicadores m_ij
for i=2:n+1
    for j=1:i-1
        BDA(i,j)= -(n+2-i)/(1-t);
    end
end 



% % %Cálculo de los multiplicadores tilde m^ij

for i=1:n+1

    for j=i+1:n+1
        BDA(i,j)= (nchoosek(n,n+1-j)*t)/(nchoosek(n,n+2-j)*(1-t));
    end
end
% 

% %Cálculo de los pivots

for i=0:n
    BDA(i+1,i+1)=factorial(i)*nchoosek(n,i)*(1-t)^((n-2*i));
end 

BDA=funcBDA(n,t);


Px=transpose(TNSolve(BDA,xb));
Py=transpose(TNSolve(BDA,yb));

PMx=A\transpose(xb);
PMy=A\transpose(yb);


P= [Px Py];
PM=[PMx PMy];

l=100;
for i=1:l
    t1(i)=i/(l+1);
end 
circunferencia_real_polinomio(t1,l)

 hold on   
%     
dibujar(P,l,t1)

dibujar(PM,l,t1)






% 
% 
% 
% 
% 
% 
