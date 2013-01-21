% REPNUM Generating repeated numbers
%
%   r = REPNUM(x, cnts);
%
%       Generating a vector r which contains each element of x, each 
%       repeated by the corresponding number of times given in cnts.
%
%       For example, 
%
%           REPNUM([10 20 30], [2 1 3])
%
%       yields: [10 10 20 30 30 30]
%
%   x = REPNUM(cnts);
%
%       This is equivalent to REPNUM(1:n, cnts) with n = length(cnts),
%       but with a slightly faster implementation.
%
%   Remarks
%   -------
%       - x is a row when cnts is a row, otherwise x is a column.
%


