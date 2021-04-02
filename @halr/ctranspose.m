function Ht = ctranspose(H)
%TRANSPOSE Conjugate transpose of the H-matrix H.

Ht = H;

if is_leafnode(H)
    if H.admissible
        Ht.U = H.V;
        Ht.V = H.U;
    else
        Ht.F = H.F';
    end
    Ht.sz = H.sz(2:-1:1);
else
    Ht.A11 = H.A11';
    Ht.A22 = H.A22';
    Ht.A12 = H.A21';
    Ht.A21 = H.A12';
    Ht.sz = H.sz(2:-1:1);
end

end
