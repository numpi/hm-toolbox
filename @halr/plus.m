function H = plus(H1, H2)
%PLUS Sum of two HALR matrices.

H = halr_plus(H1, H2);

if min(halrrank(H1), halrrank(H2)) ~= 0
    H = compress(H);
end

end

