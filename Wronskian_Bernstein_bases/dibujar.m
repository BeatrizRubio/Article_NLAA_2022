function dibujar(P,l,t)

xB = zeros(1,l);
yB = zeros(1,l);
n = size(P,1); %P tiene una dimension mas que la de los pesos


for i=1:l
    for j=1:n
        xB(i) = xB(i) + P(j,1)*nchoosek(n-1,j-1)* ((t(i)))^(j-1) *(1-(t(i)))^(n-j);
        yB(i) = yB(i) + P(j,2)*nchoosek(n-1,j-1)* ((t(i)))^(j-1) *(1-(t(i)))^(n-j);
    end
end
 
plot(xB,yB,'o')
