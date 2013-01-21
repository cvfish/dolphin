% DD_ENTROPY Entropy of discrete distributions
%
%   r = DD_ENTROPY(a);
%
%       Evaluates the entropy of the probability vectors in a.
%
%       If a is a vector, it returns a scalar v.
%       If a is a non-vector matrix, it returns a row of size [1, n],
%       where r(i) is the entropy for a(:,i).
%
%   r = DD_ENTROPY(a, 1);
%
%       Evaluates column-wise entropy.
%
%   r = DD_ENTROPY(a, 2);
%
%       Evaluates row-wise entropy.
%
