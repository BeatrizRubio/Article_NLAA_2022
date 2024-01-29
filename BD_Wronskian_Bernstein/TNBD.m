% function B=TNBD(A);
%
% Computes the bidiagonal decomposition of a matrix A
% by performing Neville elimination on A.
% For debugging purposes only. Will not maintain high relative accuracy.
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% 
% Written October 19, 2004

function B=TNBD(A);

C=A;
[m,n]=size(A);
for k=1:n
   for i=m:-1:k+1
      if (A(i,k)~=0)
      A(i,k)=A(i,k)/A(i-1,k);
      for j=k+1:n
         A(i,j)=A(i,j)-A(i,k)*A(i-1,j);
      end
      end
   end
end

B=tril(A);
A=C';
for k=1:min(m,n)
   for i=n:-1:k+1
      if (A(i,k)~=0)
      A(i,k)=A(i,k)/A(i-1,k);
      for j=k+1:min(n,m)
         A(i,j)=A(i,j)-A(i,k)*A(i-1,j);
      end
      end
   end
end
B=B+triu(A',1);
