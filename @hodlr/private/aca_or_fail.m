function [U, V, nrm] = aca_or_fail(Afun, m, n, tol, maxrank, nrm, debug)
%
% Construct a low-rank representation of C = Afun([1:m], [1:n])
% by means of adaptive cross approximation with partial pivoting (heuristic)
%
% ----------------------INPUT--------------------------------------------------
%
% Afun   matrix or function handle returning entries of the matrix
% m, n	 rows and cols of the matrix to compress
% tol    stopping criterion for the Adaptive Cross Approximation
% comp   re-compression at the end (optional), default = 0
% norm   estimate of the norm of the matrix, if empty or not given the estimate is computed on the fly
% debug  enables some debugging prints (optional), default = 0
%
% ----------------------Output-------------------------------------------------
%
% U,V    factors of the low-rank approximation Afun = U * V'
%
% -----------------------------------------------------------------------------

require_norm_output = (nargout > 2);

if ~exist('debug', 'var')
    debug = 0;
end

if isempty(maxrank)
    maxrank = round(min(m,n) / 2);
end

has_norm = exist('nrm', 'var') && ~isempty(nrm);

U = zeros(m, 0);
V = zeros(n, 0);
k = 1;
taken_row = [];

sample_size = ceil(sqrt(n));

% Select the first pivot with a random sampling
[ind, new_ind, b, a] = find_pivot(Afun, U, V, taken_row, m, n, sample_size);

while k < min(m,n)
        
    if debug
        fprintf('Iteration: %d Pivot at (%d,%d): %e\n', ...
                k, ind, new_ind, b(new_ind))
    end
    
    if has_norm && abs(b(new_ind)) * sqrt(n) <= nrm * tol

        % Check before exiting
        [ind, new_ind, b, a] = find_pivot(Afun, U, V, taken_row, m, n, sample_size);
        
        if can_stop(new_ind, tol, nrm, b, a)
            return;
 		end
	else
		if k > 1
			a = Afun((1:m)', new_ind) - U * V(new_ind, :)';
        end
    end


        
    if b(new_ind) ~= 0
       
        % a = Afun((1:m)', new_ind) - U * V(new_ind, :)';
        a = a / b(new_ind);

		% Update norm estimate, unless explicitly given
		tnrm = norm(a) * norm(b);
		if ~has_norm
			nrm = tnrm;
			has_norm = true;
		else
		    nrm = max(nrm, tnrm);
		end     
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
	    if ~has_norm
		   % Check before exiting
		   [ind, new_ind, b, a] = find_pivot(Afun, U, V, taken_row, m, n, sample_size);
		    
		   if abs(b(new_ind)) == 0
	        	nrm = 0;
				return;
			end
	    else
			return;
		end
    end
    
    k = k + 1;
    
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

end

function [ind, new_ind, row, col] = find_pivot(Afun, U, V, taken_row, m, n, ns)
    kk = 5;
    ks = max(1, ns - 2*kk);
    rowind = [ 1 : kk, ...
        randsample(setdiff(kk+1:m-kk, taken_row), ks), ...
        m-kk+1 : m ];
    colind = [ 1 : kk, randsample(1:n, ks), n-kk+1:n ];

	rowind = sort(rowind);
	colind = sort(colind);

    B = Afun(rowind', colind) - U(rowind, :) * V(colind, :)';

    [~, idx] = max(abs(B(:)));
    [r, c] = ind2sub(size(B), idx);

    indr = rowind(r);
	indc = colind(c);

    row = Afun(indr, 1:n) - U(indr, :) * V';
    col = Afun(1:m, indc) - U * V(indc, :)';
    [mr, new_indr] = max(abs(row));
	[mc, new_indc] = max(abs(col));
	if mr >= mc
		ind = indr;
		new_ind = new_indr;
        col = Afun(1:m, new_ind) - U * V(new_ind, :)';
	else
		ind = new_indc;
		row = Afun(ind, 1:n) - U(ind, :) * V';
    	% [~, new_ind] = max(abs(row));
        new_ind = indc;
        % col = Afun(1:m, ind) - U * V(ind, :)';
	end
end

function b = can_stop(new_ind, tol, nrm, row, col)
    b = norm(row) * norm(col) / abs(row(new_ind)) < tol * nrm;
end
