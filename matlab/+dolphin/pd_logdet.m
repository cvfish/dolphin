% PD_LOGDET Log-determinant of a positive definite matrix
%
%   v = PD_LOGDET(A);
%
%       Evaluates and returns the log-determinant of A.
%
%       Here, A must be a positive definite matrix.
%
%   Remarks
%   -------
%   - When A is a scalar, or an 2-by-2 or 3-by-3 matrix, the 
%     log-determinant is directly evaluated using a simple formula. 
%     For larger matrices, the evaluation is based on Cholesky 
%     decomposition of A.
%
