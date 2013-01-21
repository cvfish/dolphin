% LOGSUMEXP Logarithm of sum of exponentials
%
%   r = LOGSUMEXP(a);
%
%       Evaluates the the log(sum(exp(a))) in a numerically stable way.
%
%       If a is a vector, r will be a scalar, otherwise, r will be 
%       a row of length n, such that r(i) is the result for a(:,i).
%       