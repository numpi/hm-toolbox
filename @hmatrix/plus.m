function H = plus(H1, H2)
%PLUS Sum of two HMATRIX matrices.

H = hmatrix_plus(H1, H2);
H = compress(H);

end

