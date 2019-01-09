function H = plus(H1, H2)
%PLUS Sum two H matrices.

if ~check_cluster_equality(H1, H2)
    error('H1 + H2: Cluster or dimension mismatch in H1 and H2');
end

H = hmatrix_plus(H1, H2);
H = compress_hmatrix(H);

end

