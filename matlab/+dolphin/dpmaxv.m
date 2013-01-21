function r = dpmaxv(K, x, dim, I)
%DPMAXV Label-based dispatched maximum
%
%   r = DPMAXV(K, x, dim, I);
%
%       Evaluates maximum of rows/columns in x depending on I.
%
%       Suppose the size of x is [m, n], then:
%
%       if dim == 1, then the output r is a K-by-n matrix, such that 
%
%           r(k, :) equals max(x(I==k, :), [], 1);
%
%       if dim == 2, then the output r is an m-by-K matrix, such that
%
%           r(:, k) equals max(x(:, I==k), [], 2);
%

code = int32(2);
r = dpaccumv(code, dim, K, x, int32(I)-1);
