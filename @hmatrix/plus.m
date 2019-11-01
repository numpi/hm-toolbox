function H = plus(H1, H2)
%PLUS Sum of two HMATRIX matrices.

H = hmatrix_plus(H1, H2);

if min(hmatrixrank(H1), hmatrixrank(H2)) ~= 0
    H = compress(H);
end

end

