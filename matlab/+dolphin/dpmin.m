function r = dpmin(siz, x, I, J)
%DPMIN Label-based dispatched minimum
%
%   r = DPMIN(K, x, I);
%
%       Evaluates minimum of values in x depending on the values in I.
%
%       The output r is a column of length K, such that r(k) is the
%       minimum of x(I == k).
%
%   r = DPMIN([m, n], x, I, J);
%
%       Evaluates minimum of values in x depending on the subscripts in 
%       I and J.
%
%       The output r is an m-by-n matrix, such that r(i, j) is the minimum
%       of the x(I == i & J == j).
%
%   Remarks
%   -------
%   - Entries in r that correspond to no values in x are set to inf.
%

narginchk(3, 4)
code = int32(3);
if nargin == 3
    r = dpaccum(code, siz, x, int32(I)-1);
else
    r = dpaccum(code, siz, x, int32(I)-1, int32(J)-1);
end
