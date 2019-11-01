function H = mldivide(H1, H2)
    if is_leafnode(H1) || iszero(H1.A12)
        H = solve_lower_triangular(H1, H2);
    elseif iszero(H1.A21)
        H = solve_upper_triangular(H1, H2);
    else
        [L, U] = lu(H1);        
        H = L \ (U \ H2);
    end
end
