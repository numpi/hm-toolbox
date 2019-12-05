function [U, V] = aca_or_fail(Afun, m, n, tol, maxrank)
%ACA_OR_FAIL Try to get a low-rank representation of A, or bail our if it
%does not appear to be low-rank. 

%
% Construct a low-rank representation of C = 1 / (x(i) + y(j))
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

if isempty(maxrank)
    maxrank = round(min(m,n) / 10);
end

if ~exist('debug', 'var')
	debug = 0;
end

U = zeros(m, 0);
V = zeros(n, 0);
k = 1;
ind = 1;

% Select the first pivot
for j = 1 : 50
    first_indices = unique(randi(m, 4, 1));
    rows = Afun(first_indices, 1:n);

    [mm, ind] = max(max(abs(rows), [], 2));
    
    if mm > 0
        break;
    end
end

% Store residual history
res = [];

while k < min(m,n)
	b = Afun(ind, 1:n) - U(ind, :) * V';
	[~, new_ind] = max(abs(b));
	if debug
		fprintf('Iteration: %d Pivot at (%d,%d): %e\n', k, ind, new_ind, b(new_ind))
	end
	a = Afun((1:m)', new_ind) - U * V(new_ind, :)';
	a = a/b(new_ind);
	U = [U, a];
	V = [V, b'];
	[~, tind] = max(abs(a([1:ind - 1, ind + 1:m])));  
        if tind >= ind
		ind = tind + 1;
	else
		ind =tind;
	end

	if k == 1
		nrm = norm(a) * norm(b);
	end
	k = k + 1;
    
    res = [res, norm(a) * norm(b) ];
    
	if res(end) < tol * nrm
		break
    end
    
    if k >= maxrank || k >= min(m, n) / 2
        U = [];
        V = [];
        return;
    end
end


end

