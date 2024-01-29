%Bidiagonal decomposition  of Wronskian matrix of Bernstein basis of positive degree 
%E. Mainar, J.M. Peña, B. Rubio. Theorem 4

clear all
format longE

n=24; 

syms t

f=t;
g=1-t;
t=-1;  

%Cálculo de la matriz wronskiana

A=zeros(n+1);

for i=1:n+1
	for j=1:n+1
	A(i,j)=eval( diff(nchoosek(n,j-1)*  ((f))^(j-1) *((g))^(n+1-j),i-1));
 	end 
end

     
BDA=zeros(n+1);


%BD(JWJ) base de Bernstein

%Cálculo de los multiplicadores m^ij
for i=2:n+1
    for j=1:i-1
        BDA(i,j)= (n+2-i)/(1-t);
    end
end  

 
%Cálculo de los multiplicadores tilde m^ij

for i=1:n+1
    for j=i+1:n+1
        BDA(i,j)= (nchoosek(n,n+1-j)*t)/(nchoosek(n,n+2-j)*(t-1));
    end
end


%Cálculo de los pivots 

for i=0:n
    BDA(i+1,i+1)=factorial(i)*nchoosek(n,i)*(1-t)^((n-2*i));
end


%%Resolvemos problemas de Álgebra Lineal Numérica
J=sym(zeros(n+1,n+1));


for i=1:n+1
    J(i,i)=(-1)^(i-1);
end 


%Resolución Sistema lineal
b =[17,31,77,83,27,11,96,57,70,64,29,41,46,16,74,1,2,6,7,5,1,2,6,7,5];
SolBD=transpose(TNSolve(BDA,J*transpose(b)));
SolM = A\transpose(b)
SolB=J*transpose(SolBD)



%Cálculo de la matriz inversa
IB=TNInverseExpand(BDA);
IB2=J*IB*J %Aplicamos Corolario
IM=inv(A)
 


%Cálculo de valores propios y singualares
VPB=min(TNEigenValues_prueba(BDA));
VPM=min(eig(A))
SVB=min(TNSingularValues(BDA))
SVM=min(svd(A))
 

