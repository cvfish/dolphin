function r = dpminv(K, x, dim, I)
%DPMINV Label-based dispatched minimum
%
%   r = DPMINV(K, x, dim, I);
%
%       Evaluates minimum of rows/columns in x depending on I.
%
%       Suppose the size of x is [m, n], then:
%
%       if dim == 1, then the output r is a K-by-n matrix, such that 
%
%           r(k, :) equals min(x(I==k, :), [], 1);
%
%       if dim == 2, then the output r is an m-by-K matrix, such that
%
%           r(:, k) equals min(x(:, I==k), [], 2);
%

code = int32(3);
r = dpaccumv(code, dim, K, x, int32(I)-1);
