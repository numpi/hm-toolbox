function [U, V] = aca_or_fail(Afun, m, n, tol, maxrank)
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

if ~exist('comp', 'var')
    comp = 0;
end

if ~exist('debug', 'var')
    debug = 0;
end

if isempty(maxrank)
    maxrank = round(min(m,n) / 2);
end

U = zeros(m, 0);
V = zeros(n, 0);
k = 1;
ind = 1;

% Select the first pivot
first_indices = randsample([1:m], 10);
rows = Afun(first_indices, 1:n);

[~, ind] = max(max(abs(rows), [], 2));
ind = first_indices(ind);
taken_row = ind;
nrm = 0;
while k < min(m,n)

    b = Afun(ind, 1:n) - U(ind, :) * V';
    [~, new_ind] = max(abs(b));
    if debug
        fprintf('Iteration: %d Pivot at (%d,%d): %e\n', k, ind, new_ind, b(new_ind))
    end
    
    if abs(b(new_ind)) <= nrm * tol
    	first_indices = randsample(setdiff(1:m, taken_row), min(m - length(taken_row), 10));
        rows = Afun(first_indices, 1:n) - U(first_indices, :) * V';
        [mx, ind] = max(max(abs(rows), [], 2));
        if mx <= tol * nrm
        	return;
        else
    		ind = first_indices(ind);    
            b = Afun(ind, 1:n) - U(ind, :) * V';
            [~, new_ind] = max(abs(b));
    	end
    end
    
    a = Afun((1:m)', new_ind) - U * V(new_ind, :)';
    a = a/b(new_ind);
    U = [U, a];
    V = [V, b'];
    [~, tind] = max(abs(a([1:ind - 1, ind + 1:m])));
    if tind >= ind
        ind = tind + 1;
    else

        ind = tind;
    end

    if k == 1
        nrm = norm(a) * norm(b);
    end
    k = k + 1;
    tnrm = norm(a) * norm(b);
    nrm = max(nrm, tnrm);
    if  tnrm < tol * nrm && k < min(m,n) - 1 % If the heuristic criterion detect convergence we still perform a sample on a few rows in the residual
	first_indices = randsample(setdiff([1:m], taken_row), min(m - length(taken_row), 10));
	rows = Afun(first_indices, 1:n) - U(first_indices, :) * V';
	[mx, ind] = max(max(abs(rows), [], 2));
	if mx < tol * nrm
        	break
	else
		ind = first_indices(ind);
	end
    end
    taken_row = [taken_row, ind];

    if k >= maxrank || k >= min(m, n) / 2
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
