function H = minus(H1, H2)
%MINUS Difference of two HMATRIX matrices
H = hmatrix_minus(H1, H2);

if min(hmatrixrank(H1), hmatrixrank(H2)) ~= 0
    H = compress(H);
end

end
