function eT = expm(A, method, N, nrm)
%EXPM Evaluate the matrix exponential of A.
%
% E = EXPM(A) computes the matrix exponential of the HSS matrix A. The
% method used is a Pade' approximant combined with a scaling and squaring
% technique.
%
% E = EXPM(A, METHOD) selects a different method for the evaluation of the
% matrix exponential. The available choices are:
%
%  - 'pade': The default choice, a Pade approximant with scaling and
%    squaring
%  - 'taylor': A truncated Taylor series with scaling and squaring.

if ~exist('method','var')
    method = 'pade';
end

if ~exist('N', 'var')
    N = 26;
end

if ~exist('nrm', 'var')
    nrm = norm(A);
end

n = size(A, 2);

if nrm == 0
    eT = hss('diagonal', ones(n,1));
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
    eT = hss('diagonal', ones(n,1));
    
    tempT = eT;
    for i=1:maxit
        tempT = tempT * A * (1 / i);
        eT = eT + tempT;
    end
elseif strcmp(method,'pade')
    
    c = 1 / 2;
    eTn = hss('diagonal', ones(n,1)) + c * A;
    eTd = hss('diagonal', ones(n,1)) - c * A;
    
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
    
    % eT =  eTn / eTd;
    eT = eTd \ eTn;
else
    error('Invalid parameter method in EXPM');
end

% eT = eT^(2^h);
for i = 1 : h
    eT = eT * eT;
end

