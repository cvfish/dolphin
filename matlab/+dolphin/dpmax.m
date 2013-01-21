function r = dpmax(siz, x, I, J)
%DPMAX Label-based dispatched maximum
%
%   r = DPMAX(K, x, I);
%
%       Evaluates maximum of values in x depending on the values in I.
%
%       The output r is a column of length K, such that r(k) is the
%       maximum of x(I == k).
%
%   r = DPMAX([m, n], x, I, J);
%
%       Evaluates maximum of values in x depending on the subscripts in 
%       I and J.
%
%       The output r is an m-by-n matrix, such that r(i, j) is the maximum
%       of the x(I == i & J == j).
%
%   Remarks
%   -------
%   - Entries in r that correspond to no values in x are set to -inf.
%

narginchk(3, 4)
code = int32(2);
if nargin == 3
    r = dpaccum(code, siz, x, int32(I)-1);
else
    r = dpaccum(code, siz, x, int32(I)-1, int32(J)-1);
end
