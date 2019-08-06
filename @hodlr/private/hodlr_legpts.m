function [x, w] = hodlr_legpts(N)
%HODLR_LEGPTS Obtain the Legendre points and relative Guassian weights for [0 pi].
%
% [X, W] = LEGPTS(N) obtains the Legendre points and quadrature weights. The
% values have been precomputed using Chebfun, and are hardcoded for the
% following values of N.
%
%   N = 4, 8, 16, 32, 36, 48, 64, 80, 96, 128, 192, 256
%

path = fileparts(mfilename('fullpath'));

X = dlmread([ path, sprintf('/../../data/legpts_%d.dat', N) ], '\t');

x = X(:,1);
w = X(:,2);