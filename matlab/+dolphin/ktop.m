% KTOP Find the K smallest/largest values
%
%   r = KTOP(x, k);
%
%       Finds k smallest values of x, and returns them in ascending order.
%
%       If x is a vector, r will be a vector of length k.
%       If x is a non-vector matrix of size [m n], then r will be 
%       a matrix of size [k n], where r(:,j) corresponds to x(:,j).
%
%   [r, ri] = KTOP(x, k);
%
%       Additionally returns the indices of the selected elements.
%
%   r = KTOP(x, -k);
%   [r, ri] = KTOP(x, -k);
%
%       Finds k largest values, and returns them in descending order.
%
