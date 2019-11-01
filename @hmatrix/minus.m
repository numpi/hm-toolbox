function H = minus(H1, H2)
%MINUS Difference of two HMATRIX matrices
H = hmatrix_minus(H1, H2);
H = compress(H);
end
