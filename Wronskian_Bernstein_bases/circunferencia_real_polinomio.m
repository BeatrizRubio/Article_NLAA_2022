function circunferencia_real_polinomio(t,l)

bx1 = zeros(1,l);
by1 = zeros(1,l);
for i=1:l
    bx1(i)=cos(2*pi*t(i));
    by1(i)=sin(2*pi*t(i));
end
 
plot(bx1, by1,'+');