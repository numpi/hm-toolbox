function H = minus(H1, H2)
if isa(H1, 'hss') && isa(H2, 'hodlr')
    if ~check_cluster_equality(H1, H2)
        error('H1 - H2: Cluster or dimension mismatch in H1 and H2');
    end
    H = hss2hodlr(H1) - H2;
elseif isa(H1, 'hss') && ~isa(H2, 'hss')
    H = full(H1) - H2;
elseif ~isa(H1, 'hss') && isa(H2, 'hss')
    H = H1 - full(H2);
else
    if ~check_cluster_equality(H1, H2)
        error('H1 - H2: Cluster or dimension mismatch in H1 and H2');
    end
    H = hss_sum(H1, hss_scalar_mul(-1, H2));
end
