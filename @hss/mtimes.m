function H = mtimes(H1, H2)

% If any of the arguments if hodlr, cast everything to hodlr
if isa(H1, 'hodlr') || isa(H2, 'hodlr')
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

if isa(H2,'hss')
    if isa(H1,'hss')
        [~, c] = cluster(H1); [r, ~] = cluster(H2);
        if ~check_cluster_equality(c, r)
            error('H1 * H2: Cluster or dimension mismatch in H1 and H2');
        end
        H = hss_mat_mat_mul(H1, H2);
    elseif isscalar(H1)
        H = hss_scalar_mul(H1, H2);
    else
        H = hss_vec_mat_mul(H1, H2);
    end
else
    if isscalar(H2)
        H = hss_scalar_mul(H2, H1);
    else
        H = hss_mat_vec_mul(H1, H2);
    end
end

end


