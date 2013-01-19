function A = add_diag(A, v)
%ADD_DIAG Add values to diagonal elements
%
%   B = ADD_DIAG(A, v)
%
%       Add the value(s) of v to the diagonal elements of A.
%       
%   Parameters
%   ----------
%   - A :   A square matrix.
%   - v :   The value(s) to be added to the diagonal of A.
%           v can be either a scalar or a vector.
%
%   Returns
%   -------
%   - B :   The resultant matrix.
%

%% argument checking

if ~ismatrix(A)
    error('add_diag:invalidarg', 'A must be a matrix.');
end

%% main

siz = size(A);
inds = 1 + (0:min(siz)-1) * (siz(1)+1);

if size(v, 1) > 1
    v = v.';
end
A(inds) = A(inds) + v;

