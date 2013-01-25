% LDIST L-norm distances between vectors
%
%   r = LDIST(a, b, p);
%   r = LDIST(a, b);
%
%       Computes L-p distance between corresponding vectors in a and b.
%       When p is omitted, it computes L2 distance (Euclidean).
%
%       When p is 1, 2, or inf, special implementation (much faster) is 
%       used to perform the computation.
%
%   r = LDIST(a, b, p, w);
%
%       Computes weighted L-p distance between corresponding vectors 
%       in a and b.
%
