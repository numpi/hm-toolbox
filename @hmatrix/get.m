function M = get(H, m, n)
% Equivalent to full(H(m,n))

M = hmatrix_sub(H,m,n);

end

function M = hmatrix_sub(H, mind, nind)

if is_leafnode(H)
    if H.admissible
        M = H.U(mind, :) * H.V(nind, :)';
    else
        M = H.F(mind,nind);
    end
else
    m1 = H.A11.sz(1);
    n1 = H.A11.sz(2);
    mind1 = my_intersect(m1, mind);
    nind1 = my_intersect(n1, nind);
    
    mind = mind - m1; nind = nind - n1;
    mind2 = my_intersect(H.A22.sz(1), mind);
    nind2 = my_intersect(H.A22.sz(2), nind);
    
    M = [ hmatrix_sub(H.A11, mind1, nind1) , hmatrix_sub(H.A12, mind1, nind2) ; ...
          hmatrix_sub(H.A21, mind2, nind1) , hmatrix_sub(H.A22, mind2, nind2) ];
end

end

function s = my_intersect(m, mind)
    s = mind(mind <= m);
    s = s(s > 0);
end
