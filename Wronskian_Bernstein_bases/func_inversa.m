
function inversa =func_inversa(n,x)
g = @(t) 1 - t; 
syms t;

for i=1:n+1
    for j=1:n+1  
    if j<=i
     inverseL(i,j)=(factorial(i-1)/factorial(j-1))*nchoosek(n+1-(j),i-j)/g(t)^(i-j);
    end
    end   
end
inverse1=matlabFunction(inverseL);
inversa=inverse1(x);
