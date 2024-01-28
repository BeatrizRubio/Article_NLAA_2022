format long E

clear all

%Calculamos DU=L^-1b con L^-1b caluculado en Mathematica. 

%rng(0,'twister');

n=80


t=0.1
 c=csvread('datab.csv');
 b=csvread('datainver.csv')
 
 
% inversa=func_inversa(n,t);
%  c1=inversa*c;
%  b=transpose(c1);
%  

%Calculo de la matriz DU

for i=1:n+1
  for j=i:n+1
      if i==j
         U(i,j)=1;
      else
       U(i,j)= nchoosek(n+1-(i),n+1-(j))*((t)/(1-t))^(j-i);
      end 
   end
end



for i=1:n+1
    D(i,i)=nchoosek(n,i-1)*factorial(i-1)*(1-t)^(n-2*i+2);
end


DU=D*U;

DU


% 
% %DU=csvread('data.csv');
%  %DU=matrixDU();

% % 
% % % %Calculo de la factorizacion DU
% % %  % %Cálculo de los multiplicadores m_ij
for i=2:n+1
    for j=1:i-1
        BDA(i,j)=0;
    end
end 
% 
% % % % %Cálculo de los multiplicadores tilde m^ij
% % % % % 
for i=1:n+1

    for j=i+1:n+1
        BDA(i,j)= (nchoosek(n,n+1-j)*t)/(nchoosek(n,n+2-j)*(1-t));
    end
end
% % % 
% % % %Cálculo de los pivots
% % % 
for i=0:n
    BDA(i+1,i+1)=factorial(i)*nchoosek(n,i)*(1-t)^((n-2*i));
end 
% % 
% % 
SolB=TNSolve(BDA,b)
SolM = DU\b

%  

 dlmwrite('sistemaDUB.csv',SolB,'precision','%.100f');
 dlmwrite('sistemaDUM.csv',SolM,'precision','%.100f');
 dlmwrite('datainver.csv',b,'precision','%.100f');
%  
b
