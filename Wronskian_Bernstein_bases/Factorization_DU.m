format long E
clear all


%rng(0,'twister');

n=180
t=0.01
 r=randi([1 10],1,n+1);
 
 %r=10*rand([1 n+1]);


for i=1:n+1
    b(i)=(-1)^(i+1)*r(i);
end    
   
 
 
 inversa=func_inversa(n,t);
 c1=inversa*transpose(b);
 b=transpose(c1);
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


% 
for i=1:n+1
    D(i,i)=nchoosek(n,i-1)*factorial(i-1)*(1-t)^(n-2*i+2);
end


 DU=D*U;
% 

%DU=csvread('data.csv');
 %DU=matrixDU();

% %nuevabmatlab2=vpa(nuevabmatlab);
% 
% 

% 
% % % %Calculo de la factorizacion DU
% % %  % %Cálculo de los multiplicadores m_ij
for i=2:n+1
    for j=1:i-1
        BDA(i,j)=0;
    end
end 
% % 
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
% 
% % 
SolB=transpose(TNSolve(BDA,b))

SolM = DU\transpose(b)

b=transpose(b)
% %  
% % 
 dlmwrite('sistemaDUB.csv',SolB,'precision','%.100f');
 dlmwrite('sistemaDUM.csv',SolM,'precision','%.100f');
 dlmwrite('b.csv',b,'precision','%.100f');


