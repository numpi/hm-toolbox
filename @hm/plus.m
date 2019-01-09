function H = plus(H1, H2)
%PLUS Sum two H matrices.

%hmatrix_pack_plus(H1, H2)

if ~check_cluster_equality(H1, H2)
    error('Cluster or dimension mismatch in A and B');
end

H = hmatrix_plus(H1, H2);
H = compress_hmatrix(H);

end

