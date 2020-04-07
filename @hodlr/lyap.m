function X = lyap(varargin)
%LYAP Solve the Lyapounov equation AX + XA' + C = 0
%
% X = LYAP(A, C) solves the Lyapunov equation using a divide and conquer
%     method.
%
% X = LYAP(A, C, SA) uses the D&C method exploiting the sparse version of
%     A, passed as the argument SA, to construct the Krylov spaces.
%
% X = LYAP(A, C, 'method', METHOD) uses the method specified. The available
%     options are:
%      - 'd&c': The D&C approach of [2].
%      - 'sign': The sign function iteration, see [1,3] (works only if
%         A has eigenvalue on the open left or right half-planes).
%      - 'expm': Uses a quadrature formula to approimate the solution [3].
%         Works for positive definite matrices.
%
% X = LYAP(A, C, sA) uses the D&C approach exploiting the sparse structure
%     of A, stored in the variable SA.
%
% [1] Grasedyck, L., Hackbusch, W., & Khoromskij, B. N. (2003). Solution
%     of large scale algebraic matrix Riccati equations by use of
%     hierarchical matrices. Computing, 70(2), 121-165.
%
% [2] Kressner, D., Massei, S., & Robol, L. (2017). Low-rank updates
%     and a divide-and-conquer method for linear matrix equations.
%     arXiv preprint arXiv:1712.04349.
%
% [3] Massei, S., Palitta, D., & Robol, L. (2018). Solving Rank-Structured
%     Sylvester and Lyapunov Equations. SIAM Journal on Matrix Analysis and
%     Applications, 39(4), 1564-1590.

N = 64;

p = inputParser;

A = varargin{1};
if length(varargin) == 2 || issparse(varargin{3}) || ischar(varargin{3})
    C = varargin{2};
    is_lyapunov = true;
else
    C = varargin{3};
    B = varargin{2};
    is_lyapunov = false;
end

if isa(C, 'hmatrix')
    X = dac_lyap_blr(varargin{:});
    return;
end


addParameter(p, 'debug', false, @islogical);
addParameter(p, 'expm',  'pade', @ischar);
addParameter(p, 'method', 'd&c', @ischar);
addParameter(p, 'parallel', false, @islogical);
addParameter(p, 'autosplit', false, @islogical);

for first_keyword = 1 : length(varargin)
    if ischar(varargin{first_keyword})
        first_keyword = first_keyword - 1;
        break;
    end
end

parse(p, varargin{first_keyword+1:end});

debug = p.Results.debug;
expm_method = p.Results.expm;
method = p.Results.method;
parallel = p.Results.parallel;

if ( is_lyapunov && length(varargin) >= 3 && issparse(varargin{3})  )|| ...
        ( ~is_lyapunov && length(varargin) >= 4 && issparse(varargin{4}) )
    
    if is_lyapunov
        nrmA = 1 / norm(A);
        X = sparse_dac_lyap(A * nrmA, A' * nrmA, C * nrmA, ...
            varargin{3} * nrmA, varargin{3}' * nrmA);
    else
        nrm = 1 / max(norm(A), norm(B));
        X = sparse_dac_lyap(A * nrm, B * nrm, C * nrm, ...
            varargin{4} * nrm, varargin{5} * nrm);
    end
    
    return;
end

if strcmp(method, 'd&c')
    if is_lyapunov
        nrmA = 1 / norm(A);
        X = dac_lyap(A * nrmA, A' * nrmA, C * nrmA);
    else
        nrm = 1 ./ max(norm(A), norm(B));
        X = dac_lyap(A * nrm, B * nrm, C * nrm);
    end
    
    X = compress_hodlr(X);
    return;
end

if strcmp(method, 'sign')
    if is_lyapunov
        X = sign_lyap(A, C);
    else
        X = sign_lyap(A, B, C);
    end
    
    return;
end

if strcmp(method, 'adi')
    X = adi_lyap(A, C, 1e-6, 100);
    return;
end

[x, w] = hodlr_legpts(N);

% Acceleration parameter.
L = 5;

% The computation of the norm is required only for scaling and squaring
% methods, therefore we can save some time in the ratcheb case.
if ~strcmp(expm_method, 'ratcheb')
    nrm = norm(A);
    
    if nrm > 5.37
        A = A / (nrm / 5.37);
        C = C / (nrm / 5.37);
        nrm = 5.37;
    end
else
    nrm = 0.0;
end

% Perform the change of variables and weights
ww = (2 * L * w) .* ( sin(x) ./ (1 - cos(x)).^2 );
xx = L * cot(x ./ 2).^2;

nevals = 0;

% X = partialSum(xx, ww, @f);
if ~parallel
    X = ww(1) * f(xx(1));
    for i = 2 : length(ww)
        X = X + ww(i) * f(xx(i));
    end
else
    n = size(A, 1);
    
    try
        if exist('gcp', 'file')
            pool = gcp;
            nworkers = pool.NumWorkers * 6;
        else
            nworkers = parpool('size');
        end
    catch
        nworkers = 1;
    end
    
    X = hodlr('zeros', n, n);
    
    XX = cell(1, nworkers);
    for j = 1 : ceil(N / nworkers)
        base = (j - 1) * nworkers + 1;
        
        parfor i = 1 : min(length(ww) - base, nworkers)
            if strcmp(expm_method, 'mixed')
                if abs(xx(i+base)) * nrm > 64
                    expm_l_method = 'ratcheb';
                else
                    expm_l_method = 'pade';
                end
            else
                expm_l_method = expm_method;
            end
            
            eA = expm(-xx(i+base)*A, expm_l_method, 13, nrm * abs(xx(i+base)));
            nevals = nevals + 1;
            XX{i} = -ww(i+base) * eA * C * eA';
        end
        
        for i = 1 : nworkers
            X = X + XX{i};
        end
    end
    
    
end

    function X = partialSum(xx, ww, ef)
        %PARTIALSUM Perform the sum of the function EF at some points.
        %
        % X = PARTIALSUM(XX, WW, EF) evaluates EF at the points in XX, and
        % combines the evaluations using the weights in WW.
        if length(xx) == 1
            X = ww * ef(xx);
        else
            mp = ceil(length(xx) / 2);
            
            X = partialSum(xx(1:mp), ww(1:mp), ef) + ...
                partialSum(xx(mp+1:end), ww(mp+1:end), ef);
        end
    end

    function Y = f(x)
        if strcmp(expm_method, 'mixed')
            if abs(x) * nrm > 64
                expm_l_method = 'ratcheb';
            else
                expm_l_method = 'pade';
            end
        else
            expm_l_method = expm_method;
        end
        
        %if strcmp(expm_l_method, 'ratcheb')
        %	  fprintf('r ');
        %else
        %	  fprintf('p ');
        %end
        
        eA = expm(-x*A, expm_l_method, 13, nrm * abs(x));
        nevals = nevals + 1;
        Y = -eA * C * eA';
        
        if debug
            fprintf (' :: %d out of %d evaluations performed\n', nevals, N);
        end
    end

end

