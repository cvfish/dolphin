% RAND_PRNG Initialize or set the seed of Dolphin's random stream
%
%   RAND_PRNG
%
%       Initializes Dolphin's random stream (if it has not been
%       initialized) with default seeding.
%
%   RAND_PRNG(seed)
%
%       Set a specific seed to Dolphin's random stream.
%
%
%   Remarks
%   -------
%       All Dolphin's functions that use random numbers will generate
%       random numbers from Dolphin's underlying stream, which is not
%       the same as MATLAB's default stream.
%
%       Those function will implicitly call rand_prng when they need
%       random numbers, and therefore end-users in general may not 
%       need to invoke this function explicitly, except when they
%       have to reproduce exactly the same sequence of random numbers.
%       
%       To produce same sequences of random numbers, one may set the
%       same seed using this function.
%
