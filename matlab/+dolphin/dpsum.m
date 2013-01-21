function r = dpsum(siz, x, I, J)
%DPSUM Label-based dispatched sum
%
%   r = DPSUM(K, x, I);
%
%       Evaluates sum of values in x depending on the values in I.
%
%       The output r is a column of length K, such that r(k) is the
%       sum of x(I == k).
%
%   r = DPSUM([m, n], x, I, J);
%
%       Evaluates the sum of values in x depending on the subscripts in 
%       I and J.
%
%       The output r is an m-by-n matrix, such that r(i, j) is the sum
%       of the x(I == i & J == j).
%

narginchk(3, 4)
code = int32(1);
if nargin == 3
    r = dpaccum(code, siz, x, int32(I)-1);
else
    r = dpaccum(code, siz, x, int32(I)-1, int32(J)-1);
end
