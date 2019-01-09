function H = minus(H1,H2)
%MINUS Difference of two HODLR matrices

if ~check_cluster_equality(H1, H2)
    error('H1 - H2: Cluster or dimension mismatch in H1 and H2');
end

if ~isa(H1,'hm') || ~isa(H2,'hm')
    error('Unsupported operation');
else
    H = hmatrix_minus(H1, H2);
    H = compress_hmatrix(H);
end

end
