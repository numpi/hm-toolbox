function H = plus(H1, H2)
%PLUS Sum two H matrices.

%hmatrix_pack_plus(H1, H2)

if size(H1, 1) ~= size(H2, 1) || size(H1, 2) ~= size(H2, 2)
    error('A + B: Dimension mismatch');
end

H = hmatrix_plus(H1, H2);
H = compress_hmatrix(H);

end

