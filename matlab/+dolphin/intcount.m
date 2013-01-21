% INTCOUNT Count integer labels
%
%   r = INTCOUNT(K, I);
%
%       Counts the number of occurrences of integer 1 to K in
%       the array I.
%
%       The output r is a column of length K, such that r(k) is 
%       the number of occurrences of k in I, i.e. r(k) = sum(I == k).
%
%   r = INTCOUNT(siz, I, J);
%
%       Counts the number of occurrences of subscript pairs in I and J.
%
%       Let siz be [m, n], then the output r will be a matrix of size
%       m-by-n, such that r(i, j) is the number of times the subscripts 
%       (i, j) appear in I and J, i.e. r(i, j) = sum(I == i & J == j).
%
%       I and J must be of the same size and the same value type.
%


