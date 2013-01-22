% SYMEIG Eigenvalues/Eigenvectors for symmetric matrices
%
%   evs = SYMEIG(a);
%
%       Evaluates the eigenvalues of a symmetric matrix a, and returns
%       them in descending order.
%
%       Let a be an n-by-n matrix, then evs is a n-by-1 column vector.
%
%   [evs, U] = SYMEIG(a);
%
%       Evaluates the eigenvalues and eigenvectors of a, and returns
%       them respectively via evs and U.
%
%       Let a be an n-by-n matrix, then evs is an n-by-1 vector comprised
%       of the eigenvalues, and U is an n-by-n matrix comprised of the
%       eigenvectors. In particular, U(:,j) is the eigenvector of evs(j).
%
%   evs = SYMEIG(a, op);
%   [evs, U] = SYMEIG(a, op);
%
%       Chooses selected method to evaluate the eigen-system.
%
%       op is a character, which can take either of the following values:
%
%       - 'N':      traditional method (the one that MATLAB itself uses)
%                   This is the slowest in general.
%
%       - 'D':      Divide-and-conquer approach. (Sometimes very fast)
%
%       - 'R':      Use Relatively Robust Representations
%                   (The recommended one in LAPACK's manual)
%
%       If op is omitted, the function uses 'R' as the default choice.
%
