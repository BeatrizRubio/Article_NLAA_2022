%Cálculo de la factorización de la matriz wronkiana de la base de
%Bernstein con la factorizacion LDU (resulta ser la misma que BDA con los multipliers)
format long E
clear all

n=160

syms t

 f=t
 g=1-t
 
 t=0.1

% % % % %Cálculo de la matriz wronskiana
% % % % 
 A=zeros(n+1);
% 
for i=1:n+1
	for j=1:n+1
		A(i,j)=eval( diff(nchoosek(n,j-1)*((f))^(j-1) *((g))^(n+1-j),i-1));
 	end 
end
% % % 
% % 
BDA=zeros(n+1);

% % % Cálculo de los multiplicadores m_ij
for i=2:n+1
    for j=1:i-1
        BDA(i,j)= -(n+2-i)/(1-t);
    end
end 

% %  % %Cálculo de los multiplicadores m_ij

% % % %Cálculo de los multiplicadores tilde m^ij
% % % % 
for i=1:n+1

    for j=i+1:n+1
        BDA(i,j)= (nchoosek(n,n+1-j)*t)/(nchoosek(n,n+2-j)*(1-t));
    end
end
% % 
% % % %Cálculo de los pivots
% % % 
for i=0:n
    BDA(i+1,i+1)=factorial(i)*nchoosek(n,i)*(1-t)^((n-2*i));
end 

 %Calculo de la matriz inversa
 IB=TNInverseExpand(BDA);
  dlmwrite('inversaBDU.csv',IB,'precision','%.45f');
  IM=inv(A);
  dlmwrite('inversaMDU.csv',IM,'precision','%.45f');
% 

