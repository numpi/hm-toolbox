function H = minus(H1, H2)
%MINUS Difference of two HALR matrices
H = halr_minus(H1, H2);

if min(halrrank(H1), halrrank(H2)) ~= 0
    H = compress(H);
end

end
