function f=initialize(t,k,l);
  global Lmax Lp sh n lma

  if k<=Lp
     m=Lmax(k)-1;
     if (k>1)
        m=min(m,part(l,k-1));
     end
     for i=1:m
        l=l+lma(k);
        sh(l+1,1:t,1:n)=-ones(t,n);
        initialize(t,k+1,l);
     end
  end
