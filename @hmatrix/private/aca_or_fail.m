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

U = zeros(m, 0);
V = zeros(n, 0);
k = 1;
ind = 1;
taken_row = [];

sample_size = 50;

m_min = 128;
n_min = 128;

if m > m_min && n > n_min
	% first try with full pivoting on a subsample
    chebm = -cos((0:m_min-1)./(m_min-1)*pi);
    chebn = -cos((0:n_min-1)./(n_min-1)*pi);
	isub = unique(round(1 + (chebm + 1) ./ 2 * (m-1)));
	jsub = unique(round(1 + (chebn + 1) ./ 2 * (n-1)));
    
	[U, V, I, J] = aca_full_pivoting(Afun(isub', jsub), length(isub), length(jsub), tol);
	if length(I) >= maxrank
        	U = [];
        	V = [];
        	return;
	else
		I = isub(I);
		J = jsub(J);
        
        [UU,SS,VV] = svd(Afun(I,J));
        ISS = diag(1 ./ sqrt(diag(SS)));        
		U = Afun(1:m, I) * UU * ISS;
		V = Afun(J, 1:n)' * VV * ISS;
		[~, RU] = qr(U, 0); 
    	[~, RV] = qr(V, 0);
		nrm = norm(RU * RV');
		taken_row = I;
		k = k + size(U, 2);
		% Then, continue with partial pivoting
		%[~, ind] = max(abs(U([1:I(end) - 1, I(end) + 1:m], end)));
    		%if ind >= I(end)
        	%	ind = ind + 1;
    		%else
        	%	ind = ind;
    		%end
    	end
end

% Select the first pivot with a random sampling
first_indices = randsample(setdiff(1:m, taken_row), min(m, sample_size));
rows = Afun(first_indices, 1:n) - U(first_indices, :) * V';
[~, ind] = max(max(abs(rows), [], 2));
ind = first_indices(ind);
taken_row = [taken_row, ind];

while k < min(m,n)

    b = Afun(ind, 1:n) - U(ind, :) * V';
    [~, new_ind] = max(abs(b));
    if debug
        fprintf('Iteration: %d Pivot at (%d,%d): %e\n', k, ind, new_ind, b(new_ind))
    end
    
    if exist('nrm', 'var') && abs(b(new_ind)) <= nrm * tol
    	first_indices = randsample(setdiff(1:m, taken_row), min(m - length(taken_row), min(m - length(taken_row), sample_size)));
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
    if b(new_ind) ~= 0
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
    tnrm = norm(a) * norm(b);
    nrm = max(nrm, tnrm);
    if  tnrm < tol * nrm && k < min(m,n) - 1 % If the heuristic criterion detect convergence we still perform a sample on a few rows in the residual
	first_indices = randsample(setdiff([1:m], taken_row), min(m - length(taken_row), sample_size));
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
