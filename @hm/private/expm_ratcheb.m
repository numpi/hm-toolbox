function X = expm_ratcheb(A, m)
%EXPM_RATCHEB Implementation of the matrix exponential
%
% X = EXPM_RATCHEB(A) computes the matrix exponential of an H-matrix A,
%     using a rational approximation that is accurate only for
%     negative-definite matrices.
%
% X = EXPM_RATCHEB(A, M) performs the computation using M terms. The only
%     valid choices for M at the moment are 7, 10, 14.
%
% Author: Leonardo Robol, based on code provided by Davide Palitta.

if ~exist('m', 'var')
    m = 7;
end

n = size(A, 1);

I = hm('diagonal', ones(n, 1));
X = hm('diagonal', zeros(n, 1));

A = -A;

% compute the nodes and wheights of the rational Chebyshev expansion
% REMARK: they don't depend on t
[omega,xi] = rational(m,1);
if mod(m,2)==0
    nu=m;
else
    nu=m-1;
end

for i=1:2:nu
    AInv = inv(A - hm('diagonal', xi(i) * ones(n, 1)));
    
    % the poles come complex conjugate
    X = X + real(2 * omega(i) * AInv);
end

% compute the last term if necessary
if nu~=m
    AInv = inv(A-hm('diagonal', xi(m) * ones(n,1)));
    X = X+omega(m) * AInv;
end

