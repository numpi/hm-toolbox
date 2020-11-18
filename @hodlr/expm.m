function eT = expm(A, method, N, nrm)
%EXPM Evaluate the matrix exponential of A.
%
% E = EXPM(A) computes the matrix exponential of the HODLR matrix A. The
% method used is a Pade' approximant combined with a scaling and squaring
% technique.
%
% E = EXPM(A, METHOD) selects a different method for the evaluation of the
% matrix exponential. The available choices are:
%
%  - 'pade': The default choice, a Pade approximant with scaling and
%    squaring
%  - 'taylor': A truncated Taylor series with scaling and squaring.
%  - 'ratcheb': A rational Chebyshev approximation that only works for
%    negative definite matrices.
%
% E = EXPM(A, METHOD, N) uses a degree N approximation. The meaning of the
% integer N is method-dependent, but typically gives the number of terms
% used in the approximation.
%
% E = EXPM(A, METHOD, N, NRM) gives an estimate for the norm of A that is
% used for the scaling and squaring.

if ~exist('method','var')
    method = 'pade';
end

if ~exist('N', 'var')
    N = 26;
end

% Implementation of an efficient evaluation for negative definite matrices
if strcmp(method, 'ratcheb')
    eT = expm_ratcheb(A);
    return;
end

if strcmp(method, 'fasttaylor')
    A_pwr = nrm.pwr;
    A_k = nrm.k;
    t   = nrm.t;
    % eT = expm_fasttaylor(t, A_pwr, A_k, A, nrm.nrm, N);
    error('fasttaylor has not been implemented yet');
    return;
end

if ~exist('nrm', 'var')
    nrm = norm(A);
end

n = size(A, 2);
if nrm == 0
    eT = hodlr('diagonal', ones(n,1));
    return;
end

if strcmp('method', 'pade')
    h = max(ceil(log2(nrm / 5.37)), 0);
else
    h = max(ceil(log2(nrm / 5.37)), 0);
end

A = A * (1 / 2^h);

if strcmp(method,'taylor')
    maxit = N;
    eT = hodlr('diagonal', ones(n,1), 'cluster', cluster(A));
    
    tempT = eT;
    for i=1:maxit
        tempT = tempT * A * (1 / i);
        eT = eT + tempT;
    end
elseif strcmp(method,'pade')
    
    c = 1 / 2;
    eTn = hodlr('diagonal', ones(n,1), 'cluster', cluster(A)) + c * A;
    eTd = hodlr('diagonal', ones(n,1), 'cluster', cluster(A)) - c * A;
    
    q = floor(N / 2);
    p = 1;
    X = A;
    for k = 2 : q
        c = c * (q-k+1) / (k*(2*q-k+1));
        X = A * X;
        cX = c*X;
        eTn = eTn + cX;
        if p
            eTd = eTd + cX;
        else
            eTd = eTd - cX;
        end
        p = ~p;
    end
    
    eT =  eTn / eTd;
else
    error('Invalid parameter method in EXPM');
end

% eT = eT^(2^h);
for i = 1 : h
    eT = eT * eT;
end

