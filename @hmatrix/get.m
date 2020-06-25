function M = get(H, m, n)
% Equivalent to full(H(m,n))

M = hmatrix_sub(H,m,n);

end

function M = hmatrix_sub(H, mind, nind)

if isempty(H.A11)
    if H.admissible
        M = H.U(mind, :) * H.V(nind, :)';
    else
        M = H.F(mind,nind);
    end
else
    m1 = H.A11.sz(1);
    n1 = H.A11.sz(2);
    
    [mind1, mind2] = my_intersect(m1, mind);
    [nind1, nind2] = my_intersect(n1, nind);
    
    M11 = hmatrix_sub(H.A11, mind1, nind1);
    M12 = hmatrix_sub(H.A12, mind1, nind2);
    M21 = hmatrix_sub(H.A21, mind2, nind1);
    M22 = hmatrix_sub(H.A22, mind2, nind2);
    
    M = [ M11 , M12  ; M21, M22 ];
end

end

function [mind1, mind2] = my_intersect(m, mind)
    s = mind <= m;
    mind1 = mind(s == 1);
    mind2 = mind(s == 0) - m;
end
