function [U, V, nrm] = aca_or_fail(Afun, m, n, tol, maxrank, nrm)
%
% Construct a low-rank representation of C = Afun([1:m], [1:n])
% by means of adaptive cross approximation with partial pivoting (heuristic)
%
% ----------------------INPUT--------------------------------------------------
%
% Afun   matrix or function handle returning entries of the matrix
% m, n	 rows and cols of the matrix to compress
% tol    stopping criterion for the Adaptive Cross Approximation
% comp   re-compression at the end (optional), default =0
% debug  enables some debugging prints (optional), default = 0
%
% ----------------------Output-------------------------------------------------
%
% U,V    factors of the low-rank approximation Afun = U * V'
%
% -----------------------------------------------------------------------------

require_norm_output = ~exist('nrm', 'var');

if ~exist('comp', 'var')
    comp = 0;
end

if ~exist('debug', 'var')
    debug = 0;
end

if isempty(maxrank)
    maxrank = round(min(m,n) / 2);
end

has_norm = exist('nrm', 'var');

U = zeros(m, 0);
V = zeros(n, 0);
k = 1;
taken_row = [];

sample_size = ceil(sqrt(n));

% Select the first pivot with a random sampling
[ind, new_ind, b] = find_pivot(Afun, U, V, taken_row, m, n, sample_size);

while k < min(m,n)
        
    if debug
        fprintf('Iteration: %d Pivot at (%d,%d): %e\n', ...
                k, ind, new_ind, b(new_ind))
    end
    
    if exist('nrm', 'var') && abs(b(new_ind)) <= nrm * tol
        % Check before exiting
        [ind, new_ind, b] = find_pivot(Afun, U, V, taken_row, m, n, sample_size);
        
        if can_stop(Afun, U, V, m, n, ind, new_ind, tol, nrm)
            return;
        end
    end
    
    if b(new_ind) ~= 0
        a = Afun((1:m)', new_ind) - U * V(new_ind, :)';
        a = a/b(new_ind);
        
        U = [U, a];
        V = [V, b'];
        
        taken_row = [ taken_row, ind ];
        
        [~, tind] = max(abs(a([1:ind - 1, ind + 1:m])));
        if tind >= ind
            ind = tind + 1;
        else
            ind = tind;
        end
        
        b = Afun(ind, 1:n) - U(ind, :) * V';
        [~, new_ind] = max(abs(b));
    else
        if k == 1 && ( ~exist('nrm', 'var') || nrm == 0.0 )
            nrm = 0;
        end
        
        return
    end
    
    % Here the norm estimated using the first vectors obtained, unless a
    % norm to use as threshold has been given by the user; this is
    % particularly useful when approximating a block of a larger matrix,
    % and truncation is desired with respect to the norm of the entire
    % matrix.
    if k == 1 && ( ~exist('nrm', 'var') || nrm == 0.0 )
        nrm = norm(a) * norm(b);
    end
    
    k = k + 1;
    
    % Update norm estimate, unless explicitly given
    tnrm = norm(a) * norm(b);
    if ~has_norm
        nrm = max(nrm, tnrm);
    end        
    
    if  tnrm < tol * nrm && k < min(m,n) - 1 % If the heuristic criterion detect convergence we still perform a sample on a few rows in the residual
        [ind, new_ind, b] = find_pivot(Afun, U, V, taken_row, m, n, sample_size);
        
        if can_stop(Afun, U, V, m, n, ind, new_ind, tol, nrm)
            break;
        end
    end
    
    if k >= maxrank || k >= min(m, n) / 2
        if require_norm_output
            [~, ru] = qr(U, 0); [~, rv] = qr(V, 0);
            nrm = norm(ru * rv');
        end
        
        U = [];
        V = [];
        return;
    end
end

if comp
    [QU, RU] = qr(U, 0); % Re-compression
    [QV, RV] = qr(V, 0);
    [U, S, V] = svd(RU * RV', 'econ');
    rk = sum(diag(S) > tol * S(1,1));
    U = QU * U(:,1:rk) * sqrt(S(1:rk,1:rk));
    V = QV * V(:,1:rk) * sqrt(S(1:rk,1:rk));
end

end

function [ind, new_ind, row] = find_pivot(Afun, U, V, taken_row, m, n, ns)
    kk = 5;
    rowind = [ 1 : kk, ...
        randsample(setdiff(kk+1:m-kk, taken_row), ns - 2 * kk), ...
        m-kk+1 : m ];
    colind = [ 1 : kk, randsample(1:n, ns-2*kk), n-kk+1:n ];

    B = Afun(rowind', colind) - U(rowind, :) * V(colind, :)';

    [~, idx] = max(abs(B(:)));
    [r, ~] = ind2sub(size(B), idx);

    ind = rowind(r);

    row = Afun(ind, 1:n) - U(ind, :) * V';
    [~, new_ind] = max(abs(row));
end

function b = can_stop(Afun, U, V, m, n, ind, new_ind, tol, nrm)
    row = Afun(ind, 1:n) - U(ind, :) * V';
    col = Afun(1:m, new_ind) - U * V(new_ind, :)';
    
    b = norm(row) * norm(col) / abs(row(new_ind)) < tol * nrm;
end
