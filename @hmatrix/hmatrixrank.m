function rk = hmatrixrank(H)
%HMATRIXRANK Obtain the maximum rank in the off-diagonal blocks of H.
%


if is_leafnode(H)
    if ~H.admissible
        rk = 0;
    else
        rk = size(H.U, 2);
    end
else
    rk = max([ hmatrixrank(H.A11), hmatrixrank(H.A12), hmatrixrank(H.A21), hmatrixrank(H.A22) ]);
end

end

