% LDIST_PW Pairwise L-norm distances between vectors
%
%   r = LDIST_PW(a, b, p);
%   r = LDIST_PW(a, b);
%
%       Computes pairwise L-p distance between columns in a and b.
%       When p is omitted, it computes L2 distance (Euclidean).
%
%       When p is 1, 2, or inf, special implementation (much faster) is 
%       used to perform the computation.
%
%   r = LDIST_PW(a, [], p);
%   r = LDIST_PW(a, []);
%
%       Computes pairwise L-p distance between columns in a.
%
%   r = LDIST_PW(a, b, p, w);
%   r = LDIST_PW(a, [], p, w);
%
%       Computes pairwise weighted L-p distance between column vectors 
%       in a and b (or in a).
%
