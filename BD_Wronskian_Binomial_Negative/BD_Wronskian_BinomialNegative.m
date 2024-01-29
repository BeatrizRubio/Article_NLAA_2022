%Bidiagonal decomposition  of Wronskian matrix of Binomial negative degree 
%E. Mainar, J.M. Peña, B. Rubio, 


clear all
format longE

n=24; 


syms t

f=t
g=1-t
t=-2 %La t tiene que ser negativa para que se cumpla el Teorema. 


%Cálculo de la matriz Wronskiana

A=zeros(n+1);

for i=1:n+1
	for j=1:n+1
	A(i,j)=eval( diff(nchoosek(n,j-1)*  ((f))^(j-1) *((g))^(n+2-j),i-1));
 	end 
end
    
BDA=zeros(n+1);


%BD(JWJ) base Binomial Negativa.

for i=2:n+1
    for j=1:i-1
        BDA(i,j)= (n+3-i)/(1-t);
    end
end  

% 
% 
% 
% % 
% % % % %Cálculo de los multiplicadores tilde m^ij
% % 
 for i=1:n+1
     for j=i+1:n+1
         BDA(i,j)= (nchoosek(n,n+1-j)*t)/(nchoosek(n,n+2-j)*(t-1));
     end
 end

 BDA
% % 
% % % %Cálculo de los pivots 
% 

for i=0:n
    BDA(i+1,i+1)=factorial(i)*nchoosek(n,i)*(1-t)^((n-2*i+1));
end


% 
% % 
 J=sym(zeros(n+1,n+1));
% 
% 
 for i=1:n+1
     J(i,i)=(-1)^(i-1);
 end 


%b =[17,31,77,83,27,11,96,57,70,64,29,41,46,16,74,1,2,6,7,5,1,2,6,7,5]

SolBD=transpose(TNSolve(BDA,J*transpose(b)));

   SolM = A\transpose(b);
%   
     SolB=J*transpose(SolBD);
 
   simplify(SolB-SolM);
% 
% 
%  dlmwrite('sistemaBinomialNegativaJB25.csv',SolB,'precision','%.45f');
%      dlmwrite('sistemaBinomialNegativaJM25.csv',SolM,'precision','%.45f');
% 
% 
% 
%Calculo de la matriz inversa
%   IB=TNInverseExpand(BDA)
%    IB2=J*IB*J %Lo hacemos según el corolario del artículo. Ya que J=inv(J)
% %   
%   dlmwrite('inversaBinomialNegativaJB25.csv',IB2,'precision','%.45f');
%     
%   IM=inv(A)
%        
%    
%   dlmwrite('inversaBinomialNegativaJM25.csv',IM,'precision','%.45f');
% 
%  
% 
%  
%  
%  
%  %dlmwrite('inversaM25.csv',IM,'precision','%.45f');
% 
% % % 
% % %VALORES PROPIOS Y SINGULARES
% % 
% 
% 
%  VPB=min(TNEigenValues_prueba(BDA))
%    VPM=min(eig(A))
%   dlmwrite('VPBinomialNegativaJB25.csv',VPB,'precision','%.45f');
%   dlmwrite('VPBinomialNegativaJM25.csv',VPM,'precision','%.45f');
% 
% 
% 
% 
%  SVB=min(TNSingularValues(BDA))
%  SVM=min(svd(A))
%  dlmwrite('VSBinomialNegativaJB25.csv',SVB,'precision','%.45f');
%  dlmwrite('VSBinomialNegativaJM25.csv',SVM,'precision','%.45f');
% 
% % 
% % 
