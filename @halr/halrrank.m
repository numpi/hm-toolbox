function rk = halrrank(H)
%HALRRANK Obtain the maximum rank in the off-diagonal blocks of H.
%


if is_leafnode(H)
    if ~H.admissible
        rk = 0;
    else
        rk = size(H.U, 2);
    end
else
    rk = max([ halrrank(H.A11), halrrank(H.A12), halrrank(H.A21), halrrank(H.A22) ]);
end

end

