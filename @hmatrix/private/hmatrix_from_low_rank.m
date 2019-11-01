function H = hmatrix_from_low_rank(U, V)

m = size(U, 1);
n = size(V, 1);

H = hmatrix_build_default_tree(m, n, hmatrixoption('block-size'));

H = hmatrix_low_rank_ric(H, U, V);

end

function H = hmatrix_low_rank_ric(H, U, V)

if is_leafnode(H)
    if H.admissible        
        H.U = U;
        H.V = V;
    else
        H.F = U * V';
    end
else
    [m1, n1] = size(H.A11);
    
    H.A11 = hmatrix_low_rank_ric(H.A11, U(1:m1, :), V(1:n1, :));
    H.A12 = hmatrix_low_rank_ric(H.A12, U(1:m1, :), V(n1+1:end,:));
    H.A21 = hmatrix_low_rank_ric(H.A21, U(m1+1:end,:), V(1:n1,:));
    H.A22 = hmatrix_low_rank_ric(H.A22, U(m1+1:end,:), V(n1+1:end,:));
end

end
