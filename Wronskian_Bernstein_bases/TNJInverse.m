function C=TNJInverse(B);

% function C=TNJInverse(B);
%
% if B=BD(A), computes BD(JA^{-1}J), where J=diag(-1)^i
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% September 29, 2004

n=size(B,1);

C=eye(n);

for i=1:n-1
   for j=n-i+1:n
      C=TNAddToNext(C,B(j,i+j-n),j);
   end
end

C=TNDiagonalScale(1./diag(B),C);

for i=n-1:-1:1
   for j=n:-1:n-i+1
      C=TNAddToPrevious(C',B(j+i-n,j),1,j)';
   end
end
