function H = plus(H1, H2)
%PLUS Sum of two HODLR matrices.
if isa(H1, 'hodlr') && isa(H2, 'hss')
    if ~check_cluster_equality(H1, H2)
        error('H1 + H2: Cluster or dimension mismatch in H1 and H2');
    end
    H = hodlr_plus(H1, hss2hodlr(H2));
elseif isa(H1, 'hodlr') && ~isa(H2, 'hodlr')
    H = full(H1) + H2;
elseif ~isa(H1, 'hodlr')
    H = H1 + full(H2);
else
    if ~check_cluster_equality(H1, H2)
        error('H1 + H2: Cluster or dimension mismatch in H1 and H2');
    end
    H = hodlr_plus(H1, H2);
    
end

end

