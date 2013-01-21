% NRMEXP Normalized exponentials (softmax)
%
%   r = NRMEXP(a);
%
%       Computes normalized exponentials as follows
%
%           r(i) = exp(a(i)) / sum(exp(a))
%
%       The computation is implemented in a numerical stable way to 
%       prevent overflow for large a.
%
%       When a is a vector, r is a vector of the normalized exponentials.
%       Otherwise, the computation is performed in a colwise fashion.
%
