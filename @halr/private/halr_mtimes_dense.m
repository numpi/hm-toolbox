function x = halr_mtimes_dense(H, v)
if is_leafnode(H)
    if H.admissible
        x = H.U * (H.V' * v);
    else
        x = H.F * v;
    end
else
    n1 = H.A11.sz(2);
    n2 = H.A22.sz(2);
    x = [ ...
        halr_mtimes_dense(H.A11, v(1:n1, :)) ...
        + halr_mtimes_dense(H.A12, v(n1+1:n1+n2, :)) ; ... 
        halr_mtimes_dense(H.A21, v(1:n1, :)) ... 
        + halr_mtimes_dense(H.A22, v(n1+1:n1+n2, :)) ];
end
end
