function r = dpsumv(K, x, dim, I)
%DPSUMV Label-based dispatched sum
%
%   r = DPSUMV(K, x, dim, I);
%
%       Evaluates sum of rows/columns in x depending on the values in I.
%
%       Suppose the size of x is [m, n], then:
%
%       if dim == 1, then the output r is a K-by-n matrix, such that 
%
%           r(k, :) equals sum(x(I==k, :), 1);
%
%       if dim == 2, then the output r is an m-by-K matrix, such that
%
%           r(:, k) equals sum(x(:, I==k), 2);
%

code = int32(1);
r = dpaccumv(code, dim, K, x, int32(I)-1);
