
%Bidiagonal decomposition  of Wronskian matrix of Bernstein basis of negative degree 
%E. Mainar, J.M. Peña, B. Rubio, 
clear all
format longE

n=24; 


syms t

f=t;
g=1-t;
t=1/7; %La t tiene que estar entre 0 y 1


% Cálculo de la matriz wronskiana
A=zeros(n+1);

for i=1:n+1
	for j=1:n+1
	A(i,j)=eval( diff(nchoosek(n+j-2,j-1)*  ((-f))^(j-1) *((g))^(-n+1-j),i-1));
 	end 
end

   
BDA=zeros(n+1);


%BD(WJ) Base de Bersntein de Grado negativo 

%Hallamos los multiplicadores m_{i,j}
for i=2:n+1
    for j=1:i-1
        BDA(i,j)= (n+i-2)/(1-t);
    end
end  


%Cálculo de los multiplicadores tilde m^ij. Aplicamos weighted.
%Dividimos (n+i-2,i-1)/(n+i-3,i-2). Cambiamos a WJ
 
for i=1:n+1
    for j=i+1:n+1
         BDA(i,j)= (nchoosek(n+j-2,j-1)*t)/(nchoosek(n+j-3,j-2)*(1-t));
    end
end

% % %Cálculo de los pivots 
% 
for i=1:n+1
    BDA(i,i)=(factorial(n+i-2)/factorial(n-1))*(1-t)^((-n+2-2*i));
end
                                        
J=sym(zeros(n+1,n+1));

 for i=1:n+1
     J(i,i)=(-1)^(i-1);
 end 


%Resolver sistema lineal

SolBD=transpose(TNSolve(BDA,transpose(b)));
SolM = A\transpose(b);
SolB=J*transpose(SolBD);

dlmwrite('sistemaBernsteinGradoNegativeJB25.csv',SolB,'precision','%.45f');
dlmwrite('sistemaBernsteinGradoNegativeJM25.csv',SolM,'precision','%.45f');

%%%Calculo de la matriz inversa
 IB=TNInverseExpand(BDA);
 IB2=J*IB; %Lo hacemos según el corolario del artículo. Ya que J=inv(J)
 IM=inv(A);
 
 
dlmwrite('inversaBernsteinGradoNegativeJB10.csv',IB2,'precision','%.45f');
dlmwrite('inversaBernsteinGradoNegativeJM10.csv',IM,'precision','%.45f');

% % % %VALORES SINGULARES

 SVB=min(TNSingularValues(BDA))
 SVM=min(svd(A))
 dlmwrite('VSBernsteinGradoNegativeJB15.csv',SVB,'precision','%.45f');
 dlmwrite('VSBernsteinGradoNegativeJM15.csv',SVM,'precision','%.45f');
