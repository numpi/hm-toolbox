function H = mtimes(H1, H2)
%MTIMES Matrix multiplication

% If any of the arguments if hss, cast everything to hodlr
if isa(H1, 'hss') || isa(H2, 'hss')
    % If the clusters do not match, throw an error
    [~, c] = cluster(H1); [r, ~] = cluster(H2);
    if ~check_cluster_equality(c, r)
        error('H1 * H2: Cluster or dimension mismatch in H1 and H2');
    end
    
    if isa(H1, 'hss'); H1 = hss2hodlr(H1); end
    if isa(H2, 'hss'); H2 = hss2hodlr(H2); end
    
    H = H1 * H2;
    
    return;
end

if isa(H2, 'halr')
    H = hodlr2halr(H1) * H2;
    return;
end

% Multiplication H * v
if isfloat(H2)
    if isscalar(H2) || all(size(H2) == 1)
        if H2 == 0
            [m, n] = size(H1);
            [r, c] = cluster(H1);
            H = hodlr('zeros', m, n, 'cluster', r, c);
            return;
        end
        
        if is_leafnode(H1)
            H = hodlr();
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
        H = hodlr_mtimes_dense(H1, H2);
    end
    
    return;
end

% Multiplication w' * H
if isfloat(H1)
    if isscalar(H1)
        H = H2 * H1;
    else
        H = dense_mtimes_hodlr(H1, H2);
    end
    
    return;
end

% Multiplication of two matrices
if isa(H1, 'hodlr') && ~isa(H2, 'hodlr')
    H = full(H1) * H2;
elseif ~isa(H1, 'hodlr') && isa(H2, 'hodlr')
    H = H1 * full(H2);
else
    [~, c] = cluster(H1); [r, ~] = cluster(H2);
    if ~check_cluster_equality(c, r)
        error('H1 * H2: Cluster or dimension mismatch in H1 and H2');
    end
    H = hodlr_mtimes(H1, H2);
end

