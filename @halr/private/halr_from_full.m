function H = halr_from_full(H, A)
[m, n] = size(A);

if isempty(H)
    H = halr_build_default_tree(m, n, halroption('block-size'));
end

H = halr_from_full_rec(H, A);
end

function H = halr_from_full_rec(H, A)
if is_leafnode(H)
    if H.admissible
        [H.U, H.V] = compress_matrix(A);        
    else
        H.F = full(A);
    end
else
    [m1, n1] = size(H.A11);
    [m2, n2] = size(H.A22);
    H.A11 = halr_from_full_rec(H.A11, A(1:m1, 1:n1));
    H.A12 = halr_from_full_rec(H.A12, A(1:m1, n1+1:n1+n2));
    H.A21 = halr_from_full_rec(H.A21, A(m1+1:m1+m2, 1:n1));
    H.A22 = halr_from_full_rec(H.A22, A(m1+1:m1+m2, n1+1:n1+n2));
end
end
