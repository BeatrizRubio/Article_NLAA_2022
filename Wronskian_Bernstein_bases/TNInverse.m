function C=TNInverse(B)

% function C=TNInverse(B)
%
% Computes BD(A^{-1}) given B=BD(A)
%
% Copyright (c) 2005 Plamen Koev. See COPYRIGHT.TXT for more details.
% April 28, 2005.

n=size(B,1);

C=eye(n);

for i=1:n-1
   for j=n-i+1:n
      C=TNAddToNext(C,-B(j,i+j-n),j);
   end
end

C=TNDiagonalScale(1./diag(B),C);

for i=n-1:-1:1
   for j=n:-1:n-i+1
      C=TNAddToPrevious(C',-B(j+i-n,j),1,j)';
   end
end
