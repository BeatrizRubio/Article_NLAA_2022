function BDA=funcBDA(n,t)

BDA=zeros(n+1);
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
end 