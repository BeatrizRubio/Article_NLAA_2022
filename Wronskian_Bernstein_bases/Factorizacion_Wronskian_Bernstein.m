%Cálculo de la factorización de la matriz wronkiana de la base de
%Bernstein con la factorizacion LDU (resulta ser la misma que BDA con los multipliers)
format long E
clear all

%b=[1,-2,3,-4,5,-1,2,-1,2,-5,-2,1,-3];
%b=[0.1,-0.2,0.3,-0.4,0.5,-0.1,0.2,-0.1,0.2,-0.5,-0.2,0.1,-0.3];


%b=[1,-2,3,-4,5,-1,2,-1,2,-5,-2,1,-3,-2,-1,2,-5,-2,1,-5,-2,3,-4,5,-1,2,-3,1,-2,-5,6];

n=3
%n=size(b,2)-1

syms t

 f=t
 g=1-t
 
 t=0.1
  %inversa=func_inversa(n,t);
 %nuevabmatlab=inversa*transpose(b);
 
 

%Calculo de la matriz L
for i=1:n+1
    for j=1:n+1  
    if j<=i
       L(i,j)=(factorial(i-1)/factorial(j-1))*nchoosek(n+1-(j),i-j)*((-t)^(i-j)/(1-t)^(i-j));
    end
    end   
end

L
% 
% 
% % % % % %Cálculo de la matriz wronskiana
% % % % % 
%  A=zeros(n+1);
% 
% for i=1:n+1
% 	for j=1:n+1
% 	A(i,j)=eval( diff(nchoosek(n,j-1)*((f))^(j-1) *((g))^(n+1-j),i-1));
%  	end 
% end
% % % % 
% % % 
% BDA=zeros(n+1);
% % % % Cálculo de los multiplicadores m_ij
% % % % for i=2:n+1
% % % %     for j=1:i-1
% % % %         BDA1(i,j)= -(n+2-i)/(1-t);
% % % %     end
% % % % end 
% % % 
% % % 
% % % 
% % %  % %Cálculo de los multiplicadores m_ij
% for i=2:n+1
%     for j=1:i-1
%         BDA(i,j)=0;
%     end
% end 
% % % % % 
% % % % % % 
% % % % %Cálculo de los multiplicadores tilde m^ij
% % % % % 
% for i=1:n+1
% 
%     for j=i+1:n+1
%         BDA(i,j)= (nchoosek(n,n+1-j)*t)/(nchoosek(n,n+2-j)*(1-t));
%     end
% end
% % % 
% % % % %Cálculo de los pivots
% % % % 
% for i=0:n
%     BDA(i+1,i+1)=factorial(i)*nchoosek(n,i)*(1-t)^((n-2*i));
% end 
% 
%  SolBDMathematica=transpose(TNSolve(BDA,nuevabmathematica))
%  SolBDMatlab=(TNSolve(BDA,nuevabmatlab))
%   
%  SolM = A\transpose(b)
