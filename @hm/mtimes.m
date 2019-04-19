function H = mtimes(H1, H2)
%MTIMES Matrix multiplication

% If any of the arguments if hss, cast everything to hm
if isa(H1, 'hss') || isa(H2, 'hss')
	% If the clusters do not match, throw an error
	[~, c] = cluster(H1); [r, ~] = cluster(H2);
	if ~check_cluster_equality(c, r)
		error('H1 * H2: Cluster or dimension mismatch in H1 and H2');
	end
	
	if isa(H1, 'hss'); H1 = hss2hm(H1); end
	if isa(H2, 'hss'); H2 = hss2hm(H2); end
	
	H = H1 * H2;
	
	return;
end

% Multiplication H * v
if isfloat(H2)
    if isscalar(H2)
        if H2 == 0
            H = hm('diagonal', zeros(size(H1,1), 1));
            return;
        end
        
        if is_leafnode(H1)
            H = hm();
            H.F = H1.F * H2;
            H.sz = H1.sz;
        else
            H = H1;
            H.A11 = H1.A11 * H2;
            H.A22 = H1.A22 * H2;
            H.U21 = H1.U21 * H2;
            H.U12 = H1.U12 * H2;
        end
    else
        H = hmatrix_mtimes_dense(H1, H2);
    end
    
    return;
end

% Multiplication w' * H
if isfloat(H1)
    if isscalar(H1)
        H = H2 * H1;
    else
        H = dense_mtimes_hmatrix(H1, H2);
    end
    
    return;
end

% Multiplication of two matrices
if isa(H1, 'hm') && ~isa(H2, 'hm')
	H = full(H1) * H2;
elseif ~isa(H1, 'hm') && isa(H2, 'hm')
	H = H1 * full(H2);
else
	[~, c] = cluster(H1); [r, ~] = cluster(H2);
	if ~check_cluster_equality(c, r)
    		error('H1 * H2: Cluster or dimension mismatch in H1 and H2');	
	end
	H = hmatrix_mtimes(H1, H2);
end

