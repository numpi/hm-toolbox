function v = hmatrix_mtimes_dense(H1, v)
%HMATRIX_MTIMES_DENSE An H matrix times a dense matrix.

if is_leafnode(H1)
    v = H1.F * v;
else
    mp = H1.A11.sz(2);
    
    v1 = v(1:mp,:);
    v2 = v(mp+1:end,:);
    
    w1 = hmatrix_mtimes_dense(H1.A11, v1);
    w2 = hmatrix_mtimes_dense(H1.A22, v2);
    w2 = w2 + H1.U21 * (H1.V21' * v1);
    w1 = w1 + H1.U12 * (H1.V12' * v2);
    
    v = [ w1 ; w2 ];
end

end

